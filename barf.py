import yaml


with open('config.yml', 'r') as c:
    config = yaml.load(c, Loader=yaml.FullLoader)

arg = config['CLUMP_PARAMETERS']['clump_dist_params_dict']

print(type(arg), arg, arg['15240000'], type(arg['15240000']), arg['15240000'][0], type(arg['15240000'][0]))