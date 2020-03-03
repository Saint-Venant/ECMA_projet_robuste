import os
import shutil

maxTime = 3*60

data = os.getcwd() + '/../data_small/'
output = 'simulation3_withoutAntisymetry/'
instances = os.listdir(data)
commandFormat = './main -instanceName {} -maxTime {} -outputFileName {}'

instances = ['10_ulysses_3.tsp', '10_ulysses_6.tsp', '10_ulysses_9.tsp']

for instance in instances:
    print(instance)
    instanceName = data + instance
    try:
        outputFileName = output + instance[:-4] + '_results.dat'
        print(outputFileName)
        command = commandFormat.format(instanceName, maxTime, outputFileName)
        os.system(command)
    except Exception:
        print('PROBLEM')
        pass
