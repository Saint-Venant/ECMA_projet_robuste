import os
import shutil

maxTime = 3*60

data = os.getcwd() + '/../data_small/'
output = 'simulation3/'
instances = os.listdir(data)
commandFormat = './main -instanceName {} -maxTime {}'

for instance in instances:
    print(instance)
    instanceName = data + instance
    try:
        command = commandFormat.format(instanceName, maxTime)
        os.system(command)
        shutil.move('defaultSave.txt', 'simulation3/' + instance[:-4] + '_results.dat')
    except Exception:
        pass
