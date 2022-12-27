# AUTOGENERATED! DO NOT EDIT! File to edit: nbs/00_settings.ipynb (unless otherwise specified).

__all__ = ['print_settings', 'load_settings', 'load_settings_as_template', 'save_settings']

# Cell
import yaml

def print_settings(settings: dict):
    """Print a yaml settings file

    Args:
        settings (dict): A yaml dictionary.
    """
    print(yaml.dump(settings, default_flow_style=False))


def load_settings(path: str):
    """Load a yaml settings file.

    Args:
        path (str): Path to the settings file.
    """
    with open(path, "r") as settings_file:
        SETTINGS_LOADED = yaml.load(settings_file, Loader=yaml.FullLoader)
        return SETTINGS_LOADED


def load_settings_as_template(path: str):
    """Loads settings but removes fields that contain summary information.

    Args:
        path (str): Path to the settings file.
    """
    settings = load_settings(path)

    for _ in ['summary','failed']:
        if _ in settings:
            settings.pop(_)

    _ = 'prec_tol_calibrated'
    if 'search' in settings:
        if _ in settings['search']:
            settings['search'].pop(_)

    return settings


def save_settings(settings: dict, path: str):
    """Save settings file to path.

    Args:
        settings (dict): A yaml dictionary.
        path (str): Path to the settings file.
    """
    with open(path, "w") as file:
        yaml.dump(settings, file, sort_keys=False)