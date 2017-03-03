import sys

models = []

for arg in sys.argv:

    if arg != sys.argv[0]: models.append(arg)

for model in models:

    mod = open(model, 'r')

    mod_red = open(model + '_red', 'w')

    i = 1

    for line in mod:

        if i % 3 == 1: mod_red.write(line)

        i = i + 1

    mod.close(); mod_red.close()
