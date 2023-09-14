import numpy as np

def reduce_estructuras(EST, nvar, i):
    """
    Extracts the variable nvar at time step i from a structure vector.

    Parameters:
    EST: list
        List of structures containing the variable nvar.
    nvar: str
        Name of the variable to extract from the structures.
    i: int or list
        Time step index or list of time step indices to extract.

    Returns:
    z: list
        List of extracted values from the structures.
    """
    if isinstance(i, int):
        z = [EST[j][nvar][i] for j in range(len(EST))]
    elif isinstance(i, list) and len(i) == len(EST):
        z = [EST[j][nvar][i[j]] for j in range(len(EST))]
    else:
        raise ValueError("Invalid input for 'i'. It should be either an integer or a list of integers.")

    return z