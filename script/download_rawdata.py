import subprocess
import os

subprocess.run('wget "https://waseda.box.com/shared/static/wkn8x878xwxtwh7lyshr12rn20o29mvd.zip" -O ../data.zip', shell=True)
subprocess.run('unzip -d ../ -o ../data.zip', shell=True)

print('finish downloading raw data')