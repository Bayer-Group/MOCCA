from typing import List, Callable, Sequence, Dict
from numpy.typing import NDArray

import numpy as np

from mocca2.classes.component import Component
from mocca2.classes.compound import Compound


def cluster_components(
    components: List[Component],
    are_same: Callable[[Component, Component], bool],
    weights: Callable[[Component], float] | NDArray | Sequence[float] | None = None,
) -> Dict[int, Compound]:
    """
    Clusters similar components together and averages them to create compounds. The components are assigned compound IDs.

    Parameters
    ----------

    components: List[Component]
        Individual peak components

    are_same: Callable[[Component, Component], bool]
        Function that decides whether two components can be assigned to the same compound

    weights: Callable[[Component], float] | NDArray | Sequence[float] | None = None
        If weights are provided, they will be used to weight the individual components when creating compounds

    Returns
    -------

    List[Compound]
        Compounds created from averages of clustered components

    """

    N = len(components)

    # create adjacency matrix
    adjacency = np.ones([N, N]) == 0

    for i in range(N):
        for j in range(i):
            adj = are_same(components[i], components[j])
            adjacency[i, j] = adj
            adjacency[j, i] = adj

    # Get IDs of the individual subgraphs
    classes = _assign_classes(adjacency)

    # Assign IDs to the components
    for component, cls in zip(components, classes):
        component.compound_id = cls

    # Initialize weights
    compound_weights = None
    if weights is None:
        compound_weights = np.ones(N)
    elif callable(weights):
        compound_weights = np.array([weights(comp) for comp in components])
    else:
        compound_weights = np.array(compound_weights)
        assert (
            len(compound_weights) == N
        ), "The number of weights does not match the number of components"

    # Create averaged compounds
    compounds = {}
    for cls in set(classes):
        filter = [comp.compound_id == cls for comp in components]
        class_comps = [c for c, f in zip(components, filter) if f]
        class_weights = compound_weights[filter]

        elution_times = np.array([c.elution_time for c in class_comps])
        spectra = np.array([c.spectrum for c in class_comps])

        elution_time = int(
            np.sum(elution_times * class_weights) / np.sum(class_weights)
        )
        spectrum = np.sum(spectra.T * class_weights, axis=1) / np.sum(class_weights)

        compounds[cls] = Compound(elution_time, spectrum)

    return compounds


def _assign_classes(adjacency: NDArray) -> List[int]:
    """Assigns unique ID to each separate subgraph defined in the adjacency matrix. Returns the IDs"""

    ids = [-1] * adjacency.shape[0]
    next_id = 0
    neighbouring = []

    for new_idx, id in enumerate(ids):
        # find component with unassigned id
        if id != -1:
            continue

        # DFS
        ids[new_idx] = next_id
        neighbouring = [
            idx for idx, adj in enumerate(adjacency[new_idx]) if adj and ids[idx] == -1
        ]
        for n in neighbouring:
            ids[n] = next_id

        while len(neighbouring) != 0:
            curr_idx = neighbouring[-1]
            neighbouring.pop()

            new_neighbouring = [
                idx
                for idx, adj in enumerate(adjacency[curr_idx])
                if adj and ids[idx] == -1
            ]
            for n in new_neighbouring:
                ids[n] = next_id

            neighbouring += new_neighbouring

        # increment next id
        next_id += 1

    return ids
