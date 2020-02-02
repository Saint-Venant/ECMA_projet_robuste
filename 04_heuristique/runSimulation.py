import os
from localSearch import heuristique

data = os.getcwd() + '/../modified_data/'
output = 'simulation4/'
instances = os.listdir(data)
dtSearch = 45
nbSearch = 3
solving_time = dtSearch*nbSearch

for instance in instances:
    instanceName = data + instance
    outputFileName = instance[:-4] + '_results.dat'
    outputDirName = output + outputFileName

    if outputFileName not in os.listdir(output):
        value = heuristique(instanceName, dtSearch=dtSearch, nbSearch=nbSearch)
        print(value)
        with open(outputDirName, 'w') as f:
            f.write(instance + '\n')
            f.write('solving time = {}\n'.format(solving_time))
            f.write('best solution = {}\n'.format(value))
