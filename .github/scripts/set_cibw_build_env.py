# set_cibw_build_env.py
import os

version = os.environ.get('PYTHON_VERSION').replace('.', '')
env_file = os.environ['GITHUB_ENV']

with open(env_file, 'a') as f:
    f.write(f'CIBW_BUILD=cp{version}-*{os.linesep}')