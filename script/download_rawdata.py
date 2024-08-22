import subprocess
import os

subprocess.run('wget "https://waseda.box.com/shared/static/fz2r7j4igktjdw1hf2oh24cs7ant85aa.zip" -O ../data.zip', shell=True)
subprocess.run('unzip -d ../ -o ../data.zip', shell=True)

print('finish downloading raw data')