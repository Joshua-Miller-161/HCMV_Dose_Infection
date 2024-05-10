import numpy as np


def average_dicts(my_dict, my_dict2):
    my_dict_avg = {}
    for key in my_dict:
        if key in my_dict2:
            # Find the longer list
            longest_list = my_dict[key] if len(my_dict[key]) > len(my_dict2[key]) else my_dict2[key]
            # Compute the average element-wise
            avg_list = []
            for i in range(len(longest_list)):
                val1 = my_dict[key][i] if i < len(my_dict[key]) else 0
                val2 = my_dict2[key][i] if i < len(my_dict2[key]) else 0
                avg_list.append((val1 + val2) / 2)
            my_dict_avg[key] = avg_list
        else:
            print(f"Key '{key}' not found in my_dict2.")
    return my_dict_avg

# Example usage:
my_dict = {'a': [4, 6, 2]}
my_dict2 = {'a': [10, 2, 4, 1]}
my_dict_avg = average_dicts(my_dict, my_dict2)
print(my_dict_avg)  # Output will be {'a': [7, 4, 3, 0.5]}
