import os

# parameters
maxIter = 10000
maxTime = 5*60

data = os.getcwd() + '/../data/'
output = 'simulation1/'
instances = os.listdir(data)
commandFormat = './main -instanceName {} -maxIter {}' + \
                ' -maxTime {} -outputFileName {}'

for instance in instances:
    instanceName = data + instance
    outputFileName = output + instance[:-4] + '_results.txt'
    command = commandFormat.format(instanceName, maxIter, maxTime, \
                                   outputFileName)
    os.system(command)
