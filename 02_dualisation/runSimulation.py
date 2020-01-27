import os

data = os.getcwd() + '/../modified_data/'
model = 'robuste2Mod.mod'
tempOutput = 'output.dat'
output = 'simulation2/'
instances = os.listdir(data)
commandFormat = 'oplrun.exe {} {}'

for instance in instances:
    instanceName = data + instance
    #outputFileName = output + instance[:-4] + '_results.dat'
    command = commandFormat.format(model, instanceName)
    os.system(command)
    break
