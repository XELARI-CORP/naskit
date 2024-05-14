from setuptools import Extension
from setuptools.command.build_ext import build_ext
from setuptools.dist import Distribution
from pathlib import Path
import shutil
import os



ext_modules = [
    Extension('naskit.algo.levenshtein._levenshtein', 
              ['naskit/algo/levenshtein/levenshtein.c']),
]

def build():
    cmd = build_ext(Distribution({'ext_modules': ext_modules}))
    cmd.ensure_finalized()
    cmd.run()

    for output in cmd.get_outputs():
        relative_extension = os.path.relpath(output, cmd.build_lib)
        shutil.copyfile(output, relative_extension)


if __name__ == "__main__":
    build()