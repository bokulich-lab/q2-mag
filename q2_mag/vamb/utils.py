# ----------------------------------------------------------------------------
# Copyright (c) 2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
def _process_vamb_arg(arg_key, arg_val):
    """Creates a list with argument and its value to be consumed by VAMB.

    Argument names will be converted to command line parameters by
    appending a '--' prefix and concatenating words separated by a '_',
    e.g.: 'some_parameter_x' -> '--someParameterX'.

    Args:
        arg_key (str): Argument name.
        arg_val: Argument value.

    Returns:
        [converted_arg, arg_value]: List containing a prepared command line
            parameter and, optionally, its value.
    """
    if arg_key == "min_contig_len":
        arg_key_flag = "-m"
    elif arg_key == "threads":
        arg_key_flag = "-p"
    else:
        arg_key_flag = f"--{arg_key}"

    if isinstance(arg_val, bool) and arg_val:
        return [arg_key_flag]
    else:
        return [arg_key_flag, str(arg_val)]
