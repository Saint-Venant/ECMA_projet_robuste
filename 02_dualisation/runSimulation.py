import os
import shutil

data = os.getcwd() + '/../modified_data/'
model = 'robuste2Mod.mod'
tempOutput = 'output.dat'
output = 'simulation2/'
instances = os.listdir(data)
commandFormat = 'oplrun.exe {} {}'

for instance in instances:
    print(instance)
    instanceName = data + instance
    #outputFileName = output + instance[:-4] + '_results.dat'
    command = commandFormat.format(model, instanceName)
    os.system(command)
    shutil.move('output.dat', 'simulation2/' + instance[:-4] + '_results.dat')
    break
