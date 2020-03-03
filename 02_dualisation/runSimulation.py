import os
import shutil

data = os.getcwd() + '/../modified_data/'
model = 'robuste2Mod.mod'
tempOutput = 'output.dat'
output = 'simulation2/'
instances = os.listdir(data)
commandFormat = 'oplrun.exe {} {}'

for instance in instances:
    instanceName = data + instance
    outputFileName = instance[:-4] + '_results.dat'
    outputDirName = output + outputFileName

    if outputFileName not in os.listdir(output):
        print(instance)
        print()
        
        try:
            command = commandFormat.format(model, instanceName)
            os.system(command)
            shutil.move('output.dat', outputDirName)
        except Exception:
            pass
