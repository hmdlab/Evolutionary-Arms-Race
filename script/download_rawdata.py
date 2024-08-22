import subprocess
import os

#subprocess.run('curl "https://drive.usercontent.google.com/download?id=1Ww9JjF6nrhPStu8M_VILuN_r1VcS7o9O&confirm=xxx" -o ../data.zip', shell=True)
subprocess.run('wget "https://waseda.box.com/shared/static/fz2r7j4igktjdw1hf2oh24cs7ant85aa.zip" -O ../data.zip', shell=True)
subprocess.run('unzip -d ../ -o ../data.zip', shell=True)

print('finish downloading raw data')