import argparse
import datetime
import json
import os
import pprint
import sys
import logging
import types
from pathlib import Path
from argparse import Namespace as NamespaceOld


class Namespace(NamespaceOld):
    def __getitem__(self, key):
        if key in self.__dict__:
            return self.__dict__[key]
        else:
            return None

    def __setitem__(self, key, value):
        setattr(self, key, value)

    def get(self, key, default=None):
        if key in self.__dict__:
            return self.__dict__[key]
        else:
            return default

    def update(self, other):
        for k, v in other.items():
            setattr(self, k, v)
    
    def dict(self):
        return self.__dict__


class Arguments(argparse.ArgumentParser):
    """
    Compare to standard ArgumentParser, this class has following features:
    1. Can read json file as parameter.
    2. Automatically create output path if not existed.
    3. Save parameter to log file.
    """

    def __init__(self, *args, **kwargs):
        super(Arguments, self).__init__(add_help=True, *args, **kwargs)
        self.formatter_class = lambda prog: argparse.RawTextHelpFormatter(
            prog, max_help_position=100, width=200)
        self.add_argument('-parameter_file', type=argparse.FileType("r"), default=None,
                          help="Read parameter from a file in json format.")
        self.add_argument('-output_parameter', type=int, default=0)
        self.add_argument('-threads', type=int, default=1)
        self.add_argument('-debug', type=int, default=0)
        # self.add_argument('-path_output', type=str)

    def add_argument_from_dictionary(self, arguments: dict) -> None:
        """"
        Over current parameter
        """

        for item in arguments:
            if "-" + item not in self._option_string_actions:
                if isinstance(arguments[item], list):
                    self.add_argument("-" + item, nargs='*', default=arguments[item],
                                      help=f"Default value: {arguments[item]}")
                else:
                    self.add_argument(
                        "-" + item, type=type(arguments[item]),
                        default=arguments[item],
                        help=f"Default value: {arguments[item]}")
            else:
                action = self._option_string_actions["-" + item]
                action.default = arguments[item]

    def parse_args(self, print_parameter=True, auto_create_path=True, args=None):
        parsed_args = super(Arguments, self).parse_args(args=args, namespace=Namespace())

        # Deal with the json file
        if parsed_args.parameter_file:
            parameter_file = json.load(parsed_args.parameter_file)
            for item in parameter_file:
                setattr(parsed_args, item, parameter_file[item])
            raise NotImplementedError("Check the code before use this function")

        # Smart fill the path_output
        if "path_output" not in parsed_args.__dict__:
            setattr(parsed_args, "path_output", None)
        if parsed_args.path_output is None:
            if "file_output" in parsed_args.__dict__:
                parsed_args.path_output = Path(parsed_args.file_output).parent
            else:
                parsed_args.path_output = Path(sys.argv[0]).parent
            parsed_args.path_output.mkdir(parents=True, exist_ok=True)

        # Convert all pathname to absolute path, and create the directory if not existed
        for para_name in parsed_args.__dict__:
            if para_name.startswith("path_") or para_name.startswith("file_"):
                value = getattr(parsed_args, para_name)
                if isinstance(value, str):
                    value = Path(value)
                elif isinstance(value, list):
                    value = [Path(item) for item in value]
                setattr(parsed_args, para_name, value)

                # Create path/file if not existed.
                if auto_create_path:
                    if isinstance(value, Path):
                        if para_name.startswith("path_"):
                            value.mkdir(parents=True, exist_ok=True)
                        else:
                            value.parent.mkdir(parents=True, exist_ok=True)
                    elif isinstance(value, list):
                        for item in value:
                            if para_name.startswith("path_"):
                                item.mkdir(parents=True, exist_ok=True)
                            else:
                                item.parent.mkdir(parents=True, exist_ok=True)

        # Save parameter to log file
        log_output = Path(parsed_args.path_output) / \
            ("parameter-" + datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S") + '.log')
        if parsed_args.output_parameter:
            logging.basicConfig(filename=log_output, level=logging.DEBUG,
                                format='%(asctime)s %(message)s')
            logging.debug(pprint.pformat(vars(parsed_args)))

        if print_parameter:
            pprint.pprint(vars(parsed_args))

        return parsed_args


class ArgumentsOld(argparse.ArgumentParser):

    def __init__(self, *args, **kwargs):
        super(Arguments, self).__init__(add_help=True, *args, **kwargs)
        self.formatter_class = argparse.RawTextHelpFormatter

        self.parameter: dict = {
            'threads': None
        }

        self.add_argument('-para', type=str)
        self.add_argument('-threads', type=int)
        self.add_argument('-debug', type=int)
        self.add_argument('-path_output', type=str)
        self.add_argument('-output_parameter', type=bool)

    def add_parameter(self, para: dict) -> None:
        """"
        Over current parameter
        """

        for item in para:
            if "-" + item not in self._option_string_actions:
                if isinstance(para[item], list):
                    self.add_argument("-" + item, nargs='*')
                else:
                    self.add_argument("-" + item, type=type(para[item]))

        self.parameter.update(para)

    def add_argument_by_example(self, para: dict):
        self.set_defaults(**para)

    def parse_args(self, args=None, namespace=None, print_parameter=True):
        # Parameter order:
        # CMD input > json file > default
        # parameters_args_cmd_input > parameters_args_cmd_input["para"] > self.parameter
        parameters_args_cmd_input = vars(super(Arguments, self).parse_args(args))
        parameters_args_cmd_input = {
            k: parameters_args_cmd_input[k] for k in parameters_args_cmd_input
            if parameters_args_cmd_input[k] is not None}

        # Parameter from json file
        if parameters_args_cmd_input.get("para", None) is not None:
            parameters_json = json.load(open(parameters_args_cmd_input["para"], "rt"))
            if "para" in parameters_json:
                parameters_json.pop("para")

            # Merge parameters
            parameters_json.update(parameters_args_cmd_input)
            self.parameter.update(parameters_json)
        else:
            # Merge parameters
            self.parameter.update(parameters_args_cmd_input)

        # Fill default parameter
        if 'threads' not in self.parameter:
            self.parameter["threads"] = 1
        if 'output_parameter' not in self.parameter:
            self.parameter["output_parameter"] = False
        if 'path_output' not in self.parameter:
            self.parameter['path_output'] = ""

        # If path_output is not defined, guess it.
        if not self.parameter['path_output']:
            if 'file_output' in self.parameter and self.parameter['file_output']:
                self.parameter['path_output'] = os.path.dirname(self.parameter['file_output'])
            elif 'path_input' in self.parameter and self.parameter['path_input']:
                self.parameter['path_output'] = os.path.dirname(self.parameter['path_input'])
            elif 'file_input' in self.parameter and self.parameter['file_input']:
                self.parameter['path_output'] = os.path.dirname(self.parameter['file_input'])

        # Convert to absolute path.
        if not os.path.isabs(self.parameter['path_output']):
            self.parameter['path_output'] = \
                os.path.join(os.getcwd(), self.parameter['path_output'])

        # Fix path name
        if self.parameter["path_output"]:
            self.parameter["path_output"] = pathlib.Path(self.parameter["path_output"])
        if self.parameter.get("path_input", ""):
            self.parameter["path_input"] = pathlib.Path(self.parameter["path_input"])

        # Create path is not existed.
        if not os.path.exists(self.parameter['path_output']):
            os.makedirs(self.parameter['path_output'])

        # Add parameter file in output path
        log_output = os.path.join(
            self.parameter["path_output"],
            "parameter-" + datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S") + '.log')
        if self.parameter["output_parameter"]:
            logging.basicConfig(filename=log_output, level=logging.DEBUG,
                                format='%(asctime)s %(message)s')
            logging.debug(pprint.pformat(self.parameter))

        if print_parameter:
            pprint.pprint(self.parameter)
        return self.parameter


if __name__ == "__main__":
    args = Arguments()
    args.add_argument("-test", type=str)

    args.add_argument_from_dictionary({
        "test3": "def"
    })
    para = args.parse_args()
    print(para)
