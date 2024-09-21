import subprocess
import os

subprocess.run('wget "https://waseda.box.com/shared/static/4nebz8s7kmpzhkeu9zutuype4rfmb7fy.zip" -O ../additional_data.zip', shell=True)
subprocess.run('unzip -d ../ -o ../additional_data.zip', shell=True)

print('finish downloading the additional data')