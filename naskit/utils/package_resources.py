import sys
import importlib.resources



def get_package_path():
    if sys.version_info.minor>=9:
        path = importlib.resources.files('naskit')
    else:
        path = next(importlib.resources.path("naskit", "").gen)
        
    return path