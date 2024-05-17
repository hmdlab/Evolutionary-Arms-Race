import subprocess

subprocess.run('curl "https://drive.usercontent.google.com/download?id=1Ww9JjF6nrhPStu8M_VILuN_r1VcS7o9O&confirm=xxx" -o ../data.zip', shell=True)
subprocess.run('unzip -d ../ -o ../data.zip', shell=True)

print('finish downloading raw data')