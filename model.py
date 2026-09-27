"""Lumped hydrological model supplied for CIVE70020 coursework.

Source: coursework-supplied Model_B.py, December 2024.
Rainfall, PET and groundwater abstraction are in mm/h; storage is in mm.
"""

import numpy as np
from numba import njit
from numpy.typing import NDArray


@njit(cache=True)
def simulate(
    parameters: tuple[float, ...] | NDArray[np.float64],
    initial_states: tuple[float, ...],
    rainfall: NDArray[np.float64],
    potential_evapotranspiration: NDArray[np.float64],
    groundwater_abstraction: NDArray[np.float64],
) -> tuple[
    NDArray[np.float64],
    NDArray[np.float64],
    NDArray[np.float64],
    NDArray[np.float64],
    NDArray[np.float64],
    NDArray[np.float64],
]:
    """Return stormflow, baseflow, upper/lower soil stores and routing stores.

    Parameter order: Zr, n1, Ks, Ksb, Zg, n2, b, c, d, T1, T2, DT.
    Initial-state order: S1, S2, V1, V2. Discharges are mm/h and stores mm.
    Forcing arrays must be non-empty, one-dimensional and equal in length.
    The supplied integrator divides each DT-hour interval into 24 substeps.
    The first output row contains initial storage and zero discharge.

    Raises:
        ValueError: If parameter, state or forcing dimensions are inconsistent.
    """

    if len(parameters) != 12:
        raise ValueError("Expected 12 model parameters.")
    if len(initial_states) != 4:
        raise ValueError("Expected four initial storage states.")
    if (
        rainfall.ndim != 1
        or potential_evapotranspiration.ndim != 1
        or groundwater_abstraction.ndim != 1
    ):
        raise ValueError("Forcing arrays must be one-dimensional.")
    if len(rainfall) == 0:
        raise ValueError("Forcing arrays must not be empty.")
    if len(rainfall) != len(potential_evapotranspiration) or len(rainfall) != len(
        groundwater_abstraction
    ):
        raise ValueError("Forcing arrays must have equal lengths.")

    (
        root_depth,
        root_porosity,
        hydraulic_conductivity,
        recharge_coefficient,
        groundwater_depth,
        groundwater_porosity,
        recharge_exponent,
        baseflow_coefficient,
        baseflow_exponent,
        surface_time_constant,
        baseflow_time_constant,
        time_step,
    ) = parameters
    (
        initial_root_storage,
        initial_groundwater_storage,
        initial_surface_routing,
        initial_baseflow_routing,
    ) = initial_states

    n = len(rainfall)

    root_capacity = root_depth * root_porosity
    groundwater_capacity = groundwater_depth * groundwater_porosity

    root_storage_series = np.zeros((n,))
    groundwater_storage_series = np.zeros((n,))
    surface_routing_series = np.zeros((n,))
    baseflow_routing_series = np.zeros((n,))
    stormflow_series = np.zeros((n,))
    baseflow_series = np.zeros((n,))

    root_storage_series[0] = initial_root_storage
    groundwater_storage_series[0] = initial_groundwater_storage
    surface_routing_series[0] = initial_surface_routing
    baseflow_routing_series[0] = initial_baseflow_routing

    root_storage = root_storage_series[0]
    groundwater_storage = groundwater_storage_series[0]
    surface_routing = surface_routing_series[0]
    baseflow_routing = baseflow_routing_series[0]

    substeps = 24
    substep_hours = time_step / substeps

    for i in range(1, n):
        rainfall_rate = rainfall[i - 1]  # current step rainfall
        pet_rate = potential_evapotranspiration[i - 1]  # current step PET
        abstraction_rate = groundwater_abstraction[i - 1]

        stormflow_volume = 0  # initialize discharge for this time step
        baseflow_volume = 0  # initialize discharge for this time step

        for _ in range(substeps):
            actual_evapotranspiration = (root_storage / root_capacity) * pet_rate

            if root_storage == root_capacity:
                infiltration = 0
            else:
                infiltration = max(
                    min(
                        rainfall_rate,
                        hydraulic_conductivity,
                        (root_capacity - root_storage) / substep_hours,
                    ),
                    0,
                )

            surface_runoff = rainfall_rate - infiltration

            if (groundwater_storage == groundwater_capacity) or (root_storage <= 0):
                recharge = 0
            else:
                recharge = recharge_coefficient * (
                    (root_storage / root_capacity) ** recharge_exponent
                )

            baseflow = baseflow_coefficient * (
                (groundwater_storage / groundwater_capacity) ** baseflow_exponent
            )

            root_storage = root_storage + substep_hours * (
                infiltration - actual_evapotranspiration - recharge
            )  # [mm]
            groundwater_storage = groundwater_storage + substep_hours * (
                recharge - baseflow - abstraction_rate
            )  # [mm]
            surface_routing = surface_routing + substep_hours * (
                surface_runoff - (surface_routing / surface_time_constant)
            )  # [mm]
            baseflow_routing = baseflow_routing + substep_hours * (
                baseflow - (baseflow_routing / baseflow_time_constant)
            )  # [mm]

            # Remove excess water if saturated
            if root_storage > root_capacity:
                surface_routing = surface_routing + (root_storage - root_capacity)  # [mm]
                root_storage = root_capacity

            if root_storage < 0 or np.isnan(root_storage):
                root_storage = 0

            if groundwater_storage > groundwater_capacity:
                baseflow_routing = baseflow_routing + (
                    groundwater_storage - groundwater_capacity
                )  # [mm]
                groundwater_storage = groundwater_capacity

            if groundwater_storage < 0 or np.isnan(groundwater_storage):
                groundwater_storage = 0

            if surface_routing < 0 or np.isnan(surface_routing):
                surface_routing = 0

            if baseflow_routing < 0 or np.isnan(baseflow_routing):
                baseflow_routing = 0

            stormflow_volume = stormflow_volume + substep_hours * (
                surface_routing / surface_time_constant
            )
            baseflow_volume = baseflow_volume + substep_hours * (
                baseflow_routing / baseflow_time_constant
            )

        root_storage_series[i] = root_storage
        groundwater_storage_series[i] = groundwater_storage
        surface_routing_series[i] = surface_routing
        baseflow_routing_series[i] = baseflow_routing
        stormflow_series[i] = stormflow_volume / time_step
        baseflow_series[i] = baseflow_volume / time_step

    return (
        stormflow_series,
        baseflow_series,
        root_storage_series,
        groundwater_storage_series,
        surface_routing_series,
        baseflow_routing_series,
    )


def least_squares(observed: NDArray[np.float64], simulated: NDArray[np.float64]) -> float:
    """Return the sum of squared residuals for paired observations and predictions.

    Arrays must be non-empty and one-dimensional, with matching shapes.
    """
    if observed.ndim != 1 or simulated.ndim != 1:
        raise ValueError("Observed and simulated values must be one-dimensional.")
    if observed.size == 0 or observed.shape != simulated.shape:
        raise ValueError("Observed and simulated values must have matching non-empty shapes.")
    return float(np.sum((observed - simulated) ** 2))


def objective_function(
    parameters: NDArray[np.float64],
    initial_states: tuple[float, ...],
    rainfall: NDArray[np.float64],
    potential_evapotranspiration: NDArray[np.float64],
    groundwater_abstraction: NDArray[np.float64],
    observed_flow: NDArray[np.float64],
) -> float:
    """Return discharge SSE over the full calibration record, including spin-up."""
    stormflow, baseflow, _, _, _, _ = simulate(
        parameters,
        initial_states,
        rainfall,
        potential_evapotranspiration,
        groundwater_abstraction,
    )
    return least_squares(observed_flow, stormflow + baseflow)
