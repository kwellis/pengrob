"""Binary Interactions

A place to store code for looking up binary interaction parameters.
Values were found using Peng Robinson Hysys and copying into excel.
Chat GPT used to convert the excel file into a python dictionary.
Python dictionary is shown below for value lookup.
"""

# binary interaction parameter dictionary
bini_dict = {
    "c1": {
        "c1": 0.0,
        "c2": 2.24e-3,
        "c3": 6.83e-3,
        "ic4": 1.31e-2,
        "nc4": 1.23e-2,
        "ic5": 1.76e-2,
        "nc5": 1.79e-2,
        "nc6": 2.35e-2,
        "nc7": 2.89e-2,
        "nc8": 3.42e-2,
        "nc9": 3.89e-2,
        "nc10": 4.36e-2,
        "n2": 3.60e-2,
        "o2": 0.0,
        "co2": 1.00e-1,
        "h2s": 8.50e-2,
        "h2o": 5.00e-1,
        "c7+": 4.36e-2,
    },
    "c2": {
        "c1": 2.24e-3,
        "c2": 0.0,
        "c3": 1.26e-3,
        "ic4": 4.57e-3,
        "nc4": 4.10e-3,
        "ic5": 7.41e-3,
        "nc5": 7.61e-3,
        "nc6": 1.14e-2,
        "nc7": 1.53e-2,
        "nc8": 1.93e-2,
        "nc9": 2.30e-2,
        "nc10": 2.67e-2,
        "n2": 5.00e-2,
        "o2": 0.0,
        "co2": 1.30e-1,
        "h2s": 8.40e-2,
        "h2o": 5.00e-1,
        "c7+": 2.67e-2,
    },
    "c3": {
        "c1": 6.83e-3,
        "c2": 1.26e-3,
        "c3": 0.0,
        "ic4": 1.04e-3,
        "nc4": 8.19e-4,
        "ic5": 2.58e-3,
        "nc5": 2.70e-3,
        "nc6": 5.14e-3,
        "nc7": 7.89e-3,
        "nc8": 1.09e-2,
        "nc9": 1.37e-2,
        "nc10": 1.66e-2,
        "n2": 8.00e-2,
        "o2": 0.0,
        "co2": 1.35e-1,
        "h2s": 7.50e-2,
        "h2o": 5.00e-1,
        "c7+": 1.66e-2,
    },
    "ic4": {
        "c1": 1.31e-2,
        "c2": 4.57e-3,
        "c3": 1.04e-3,
        "ic4": 0.0,
        "nc4": 1.34e-5,
        "ic5": 3.46e-4,
        "nc5": 3.90e-4,
        "nc6": 1.57e-3,
        "nc7": 3.22e-3,
        "nc8": 5.21e-3,
        "nc9": 7.25e-3,
        "nc10": 9.45e-3,
        "n2": 9.50e-2,
        "o2": 0.0,
        "co2": 1.30e-1,
        "h2s": 5.00e-2,
        "h2o": 5.00e-1,
        "c7+": 9.45e-3,
    },
    "nc4": {
        "c1": 1.23e-2,
        "c2": 4.10e-3,
        "c3": 8.19e-4,
        "ic4": 1.34e-5,
        "nc4": 0.0,
        "ic5": 4.95e-4,
        "nc5": 5.47e-4,
        "nc6": 1.87e-3,
        "nc7": 3.65e-3,
        "nc8": 5.75e-3,
        "nc9": 7.88e-3,
        "nc10": 1.02e-2,
        "n2": 9.00e-2,
        "o2": 0.0,
        "co2": 1.30e-1,
        "h2s": 6.00e-2,
        "h2o": 5.00e-1,
        "c7+": 1.02e-2,
    },
    "ic5": {
        "c1": 1.76e-2,
        "c2": 7.41e-3,
        "c3": 2.58e-3,
        "ic4": 3.46e-4,
        "nc4": 4.95e-4,
        "ic5": 0.0,
        "nc5": 1.25e-6,
        "nc6": 4.40e-4,
        "nc7": 1.46e-3,
        "nc8": 2.88e-3,
        "nc9": 4.45e-3,
        "nc10": 6.21e-3,
        "n2": 9.50e-2,
        "o2": 0.0,
        "co2": 1.25e-1,
        "h2s": 6.00e-2,
        "h2o": 5.00e-1,
        "c7+": 6.21e-3,
    },
    "nc5": {
        "c1": 1.79e-2,
        "c2": 7.61e-3,
        "c3": 2.70e-3,
        "ic4": 3.90e-4,
        "nc4": 5.47e-4,
        "ic5": 1.25e-6,
        "nc5": 0.0,
        "nc6": 3.93e-4,
        "nc7": 1.37e-3,
        "nc8": 2.76e-3,
        "nc9": 4.30e-3,
        "nc10": 6.03e-3,
        "n2": 1.00e-1,
        "o2": 0.0,
        "co2": 1.25e-1,
        "h2s": 6.50e-2,
        "h2o": 4.80e-1,
        "c7+": 6.03e-3,
    },
    "nc6": {
        "c1": 2.35e-2,
        "c2": 1.14e-2,
        "c3": 5.14e-3,
        "ic4": 1.57e-3,
        "nc4": 1.87e-3,
        "ic5": 4.40e-4,
        "nc5": 3.93e-4,
        "nc6": 0.0,
        "nc7": 2.97e-4,
        "nc8": 1.07e-3,
        "nc9": 2.10e-3,
        "nc10": 3.35e-3,
        "n2": 1.49e-1,
        "o2": 0.0,
        "co2": 1.25e-1,
        "h2s": 6.00e-2,
        "h2o": 5.00e-1,
        "c7+": 3.35e-3,
    },
    "nc7": {
        "c1": 2.89e-2,
        "c2": 1.53e-2,
        "c3": 7.89e-3,
        "ic4": 3.22e-3,
        "nc4": 3.65e-3,
        "ic5": 1.46e-3,
        "nc5": 1.37e-3,
        "nc6": 2.97e-4,
        "nc7": 0.0,
        "nc8": 2.41e-4,
        "nc9": 8.18e-4,
        "nc10": 1.66e-3,
        "n2": 1.44e-1,
        "o2": 0.0,
        "co2": 1.20e-1,
        "h2s": 6.00e-2,
        "h2o": 5.00e-1,
        "c7+": 1.66e-3,
    },
    "nc8": {
        "c1": 3.42e-2,
        "c2": 1.93e-2,
        "c3": 1.09e-2,
        "ic4": 5.21e-3,
        "nc4": 5.75e-3,
        "ic5": 2.88e-3,
        "nc5": 2.76e-3,
        "nc6": 1.07e-3,
        "nc7": 2.41e-4,
        "nc8": 0.0,
        "nc9": 1.71e-4,
        "nc10": 6.36e-4,
        "n2": 1.00e-1,
        "o2": 0.0,
        "co2": 1.15e-1,
        "h2s": 5.50e-2,
        "h2o": 5.00e-1,
        "c7+": 6.36e-4,
    },
    "nc9": {
        "c1": 3.89e-2,
        "c2": 2.30e-2,
        "c3": 1.37e-2,
        "ic4": 7.25e-3,
        "nc4": 7.88e-3,
        "ic5": 4.45e-3,
        "nc5": 4.30e-3,
        "nc6": 2.10e-3,
        "nc7": 8.18e-4,
        "nc8": 1.71e-4,
        "nc9": 0.0,
        "nc10": 1.48e-4,
        "n2": 1.00e-1,
        "o2": 0.0,
        "co2": 1.01e-1,
        "h2s": 5.00e-2,
        "h2o": 5.00e-1,
        "c7+": 1.48e-4,
    },
    "nc10": {
        "c1": 4.36e-2,
        "c2": 2.67e-2,
        "c3": 1.66e-2,
        "ic4": 9.45e-3,
        "nc4": 1.02e-2,
        "ic5": 6.21e-3,
        "nc5": 6.03e-3,
        "nc6": 3.35e-3,
        "nc7": 1.66e-3,
        "nc8": 6.36e-4,
        "nc9": 1.48e-4,
        "nc10": 0.0,
        "n2": 1.32e-1,
        "o2": 0.0,
        "co2": 1.18e-1,
        "h2s": 4.50e-2,
        "h2o": 5.00e-1,
        "c7+": 1.48e-4,
    },
    "n2": {
        "c1": 3.60e-2,
        "c2": 5.00e-2,
        "c3": 8.00e-2,
        "ic4": 9.50e-2,
        "nc4": 9.00e-2,
        "ic5": 9.50e-2,
        "nc5": 1.00e-1,
        "nc6": 1.49e-1,
        "nc7": 1.44e-1,
        "nc8": 1.00e-1,
        "nc9": 1.00e-1,
        "nc10": 1.32e-1,
        "n2": 0.0,
        "o2": -1.20e-2,
        "co2": -2.00e-2,
        "h2s": 1.68e-1,
        "h2o": -3.20e-1,
        "c7+": 1.32e-1,
    },
    "o2": {
        "c1": 0.0,
        "c2": 0.0,
        "c3": 0.0,
        "ic4": 0.0,
        "nc4": 0.0,
        "ic5": 0.0,
        "nc5": 0.0,
        "nc6": 0.0,
        "nc7": 0.0,
        "nc8": 0.0,
        "nc9": 0.0,
        "nc10": 0.0,
        "n2": -1.20e-2,
        "o2": 0.0,
        "co2": 9.75e-2,
        "h2s": 0.0,
        "h2o": 0.0,
        "c7+": 0.0,
    },
    "co2": {
        "c1": 1.00e-1,
        "c2": 1.30e-1,
        "c3": 1.35e-1,
        "ic4": 1.30e-1,
        "nc4": 1.30e-1,
        "ic5": 1.25e-1,
        "nc5": 1.25e-1,
        "nc6": 1.25e-1,
        "nc7": 1.20e-1,
        "nc8": 1.15e-1,
        "nc9": 1.01e-1,
        "nc10": 1.18e-1,
        "n2": -2.00e-2,
        "o2": 0.0,
        "co2": 0.0,
        "h2s": 0.0,
        "h2o": 0.0,
        "c7+": 0.0,
    },
    "h2s": {
        "c1": 8.50e-2,
        "c2": 8.40e-2,
        "c3": 7.50e-2,
        "ic4": 5.00e-2,
        "nc4": 6.00e-2,
        "ic5": 6.00e-2,
        "nc5": 6.50e-2,
        "nc6": 6.00e-2,
        "nc7": 6.00e-2,
        "nc8": 5.50e-2,
        "nc9": 5.00e-2,
        "nc10": 4.50e-2,
        "n2": 1.68e-1,
        "o2": 0.0,
        "co2": 0.0,
        "h2s": 0.0,
        "h2o": 0.0,
        "c7+": 0.0,
    },
    "h2o": {
        "c1": 5.00e-1,
        "c2": 5.00e-1,
        "c3": 5.00e-1,
        "ic4": 5.00e-1,
        "nc4": 5.00e-1,
        "ic5": 5.00e-1,
        "nc5": 4.80e-1,
        "nc6": 5.00e-1,
        "nc7": 5.00e-1,
        "nc8": 5.00e-1,
        "nc9": 5.00e-1,
        "nc10": 5.00e-1,
        "n2": -3.20e-1,
        "o2": 0.0,
        "co2": 0.0,
        "h2s": 0.0,
        "h2o": 0.0,
        "c7+": 0.0,
    },
    "c7+": {
        "c1": 4.36e-2,
        "c2": 2.67e-2,
        "c3": 1.66e-2,
        "ic4": 9.45e-3,
        "nc4": 1.02e-2,
        "ic5": 6.21e-3,
        "nc5": 6.03e-3,
        "nc6": 3.35e-3,
        "nc7": 1.66e-3,
        "nc8": 6.36e-4,
        "nc9": 1.48e-4,
        "nc10": 1.48e-4,
        "n2": 1.32e-1,
        "o2": 0.0,
        "co2": 1.18e-1,
        "h2s": 4.50e-2,
        "h2o": 5.00e-1,
        "c7+": 0.0,
    },
}
