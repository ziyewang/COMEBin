import numpy as np


def calculate_n50(lengths):
    values = np.asarray(lengths, dtype=np.int64)
    if values.ndim != 1 or values.size == 0:
        raise ValueError('N50 requires a nonempty one-dimensional length array.')
    if np.any(values < 0):
        raise ValueError('N50 lengths must be non-negative.')

    total_length = sum(int(value) for value in values)
    if total_length <= 0:
        raise ValueError('N50 requires a positive total sequence length.')
    if total_length > np.iinfo(np.int64).max:
        raise OverflowError('The total sequence length exceeds int64.')

    ordered = np.sort(values)[::-1]
    target = (total_length + 1) // 2
    index = int(np.searchsorted(np.cumsum(ordered, dtype=np.int64), target, side='left'))
    return int(ordered[index])
