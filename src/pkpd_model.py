"""Reproduce the core Berger-Rodbard subcutaneous insulin PK model.

The script simulates depot absorption and plasma-insulin profiles for four
insulin formulations, exports summary metrics, and performs a transparent
local sensitivity analysis. It is intended for educational modeling work,
not clinical decision-making.
"""

from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp


@dataclass(frozen=True)
class Formulation:
    """Parameters for the dose-dependent absorption curve."""

    s: float
    a: float
    b: float


FORMULATIONS = {
    "Regular": Formulation(s=2.0, a=0.05, b=1.7),
    "NPH": Formulation(s=2.0, a=0.18, b=4.9),
    "Lente": Formulation(s=2.4, a=0.15, b=6.2),
    "Ultralente": Formulation(s=2.5, a=0.00, b=13.0),
}

DOSE_UNITS = 24.0
ELIMINATION_RATE_PER_HOUR = 5.4
DISTRIBUTION_VOLUME_LITERS = 12.0
TIME_HOURS = np.linspace(0.0, 24.0, 1000)


def half_absorption_time(dose: float, formulation: Formulation) -> float:
    """Return T50(D) = aD + b in hours."""

    return formulation.a * dose + formulation.b


def depot_remaining_percent(
    time: np.ndarray, dose: float, formulation: Formulation
) -> np.ndarray:
    """Return the percentage of insulin remaining at the injection site."""

    t50 = half_absorption_time(dose, formulation)
    return 100.0 * t50**formulation.s / (t50**formulation.s + time**formulation.s)


def absorption_rate(time: float, dose: float, formulation: Formulation) -> float:
    """Return the rate at which insulin leaves the subcutaneous depot."""

    if time == 0.0 and formulation.s < 1.0:
        return 0.0
    t50 = half_absorption_time(dose, formulation)
    numerator = dose * formulation.s * t50**formulation.s * time ** (formulation.s - 1.0)
    denominator = (t50**formulation.s + time**formulation.s) ** 2
    return numerator / denominator


def simulate_plasma_insulin(
    time: np.ndarray,
    dose: float,
    formulation: Formulation,
    elimination_rate: float = ELIMINATION_RATE_PER_HOUR,
    distribution_volume: float = DISTRIBUTION_VOLUME_LITERS,
) -> np.ndarray:
    """Solve the one-compartment plasma-insulin balance equation."""

    def model(current_time: float, amount: np.ndarray) -> list[float]:
        entry = absorption_rate(current_time, dose, formulation)
        return [entry - elimination_rate * amount[0]]

    solution = solve_ivp(
        model,
        (float(time[0]), float(time[-1])),
        y0=[0.0],
        t_eval=time,
        rtol=1e-8,
        atol=1e-10,
    )
    if not solution.success:
        raise RuntimeError(solution.message)
    return solution.y[0] / distribution_volume


def summarize_profile(time: np.ndarray, concentration: np.ndarray) -> dict[str, float]:
    """Calculate peak, time-to-peak, and area under the curve."""

    peak_index = int(np.argmax(concentration))
    return {
        "peak_concentration": float(concentration[peak_index]),
        "time_to_peak_hours": float(time[peak_index]),
        "auc_0_24": float(np.trapezoid(concentration, time)),
    }


def local_sensitivity(
    formulation_name: str,
    formulation: Formulation,
    perturbation: float = 0.10,
) -> list[dict[str, float | str]]:
    """Estimate normalized local sensitivity of plasma-insulin AUC."""

    baseline_curve = simulate_plasma_insulin(TIME_HOURS, DOSE_UNITS, formulation)
    baseline_auc = summarize_profile(TIME_HOURS, baseline_curve)["auc_0_24"]
    rows: list[dict[str, float | str]] = []

    for parameter_name in ("s", "a", "b"):
        base_value = getattr(formulation, parameter_name)
        if base_value == 0.0:
            rows.append(
                {
                    "formulation": formulation_name,
                    "parameter": parameter_name,
                    "normalized_sensitivity": np.nan,
                }
            )
            continue

        values = formulation.__dict__.copy()
        values[parameter_name] = base_value * (1.0 + perturbation)
        increased = Formulation(**values)
        increased_auc = summarize_profile(
            TIME_HOURS,
            simulate_plasma_insulin(TIME_HOURS, DOSE_UNITS, increased),
        )["auc_0_24"]

        values[parameter_name] = base_value * (1.0 - perturbation)
        decreased = Formulation(**values)
        decreased_auc = summarize_profile(
            TIME_HOURS,
            simulate_plasma_insulin(TIME_HOURS, DOSE_UNITS, decreased),
        )["auc_0_24"]

        elasticity = (increased_auc - decreased_auc) / (
            2.0 * perturbation * baseline_auc
        )
        rows.append(
            {
                "formulation": formulation_name,
                "parameter": parameter_name,
                "normalized_sensitivity": float(elasticity),
            }
        )

    return rows


def main() -> None:
    """Run simulations and save reproducible figures and tables."""

    output_directory = Path(__file__).resolve().parents[1] / "outputs"
    output_directory.mkdir(exist_ok=True)

    figure, axes = plt.subplots(2, 1, figsize=(9, 8), sharex=True)
    summary_rows = []
    sensitivity_rows = []

    for name, formulation in FORMULATIONS.items():
        depot = depot_remaining_percent(TIME_HOURS, DOSE_UNITS, formulation)
        concentration = simulate_plasma_insulin(TIME_HOURS, DOSE_UNITS, formulation)
        normalized_concentration = 100.0 * concentration / concentration.max()

        axes[0].plot(TIME_HOURS, depot, label=name, linewidth=2)
        axes[1].plot(TIME_HOURS, normalized_concentration, label=name, linewidth=2)

        metrics = summarize_profile(TIME_HOURS, concentration)
        summary_rows.append(
            {
                "formulation": name,
                "dose_units": DOSE_UNITS,
                "t50_hours": half_absorption_time(DOSE_UNITS, formulation),
                **metrics,
            }
        )
        sensitivity_rows.extend(local_sensitivity(name, formulation))

    axes[0].set_title("Insulin Remaining at the Subcutaneous Depot")
    axes[0].set_ylabel("Remaining dose (%)")
    axes[1].set_title("Normalized Plasma-Insulin Profiles")
    axes[1].set_xlabel("Time after injection (hours)")
    axes[1].set_ylabel("Normalized concentration (%)")
    for axis in axes:
        axis.grid(alpha=0.25)
        axis.legend()

    figure.tight_layout()
    figure.savefig(output_directory / "insulin_profiles.png", dpi=200)
    plt.close(figure)

    pd.DataFrame(summary_rows).to_csv(output_directory / "summary_metrics.csv", index=False)
    pd.DataFrame(sensitivity_rows).to_csv(
        output_directory / "local_sensitivity.csv", index=False
    )

    print(pd.DataFrame(summary_rows).to_string(index=False))


if __name__ == "__main__":
    main()

