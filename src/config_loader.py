import yaml 
import pathlib
import os
import warnings

parent_dir = pathlib.Path(__file__).resolve().parent
config_file = 'config.yaml'
DEFAULT_CONFIG_FILE = os.path.join(parent_dir, config_file)

class InvalidConfigFile(Exception):
    '''
    Any sort of problem with the config file
    '''

class ConfigParamNotFound(Exception): 
    '''
    Thrown when we cannot find a key in the config
    '''

class ConflictingConfigParameters(Warning):
    '''
    Raised when something isn't quite right but we can work around it
    '''

class Config(dict): 
    
    def __init__(self, config_file: str = DEFAULT_CONFIG_FILE): 
        '''Load in the config file, check that all feilds are valid and raise any exceptions or warnings as needed

        :param str config_file: Path to the config file
            Defaults to './config.yaml'
        '''
        if config_file.split('.')[-1] not in ['yaml', 'yml']:
            raise InvalidConfigFile('Config file must be a yaml file')

        config = pathlib.Path(config_file)

        # if config.name.split('.')[-1] not in ['yaml', 'yml']:
        #     raise InvalidConfigFile('Config file must be a yaml file')

        if not config.is_file():
            raise InvalidConfigFile(f'Config file {config_file} does not exist')
        
        self.config = yaml.safe_load(open(config_file))
        self._check_config()

    def _check_config(self):
        '''Make sure certain fields are valid 
        '''
        if not pathlib.Path(self['spectra_dir']).is_dir():
            raise InvalidConfigFile(f'Spectra directory {self["spectra_dir"]} is not a vald directory')

        if not pathlib.Path(self['database_file']).is_file():
            raise InvalidConfigFile(f'Database file {self["database_file"]} is not a valid file')

        if not pathlib.Path(self['output_dir']).is_dir():
            raise InvalidConfigFile(f'Output directory {self["output_dir"]} is not a valid directory')

        if self['num_peaks'] > 0 and self['relative_abundance'] > 0: 
            warnings.warn('The two peak filtering criteria are both set to non-zero values.\n'
                                              f'Defaulting to the number of peaks set to {self["num_peaks"]}', 
                                              ConflictingConfigParameters)

    def _finditem(self, obj, key):
        '''Recursivley search for a key
        '''
        if key in obj: return obj[key]

        for _, v in obj.items():
            if isinstance(v,dict):
                item = self._finditem(v, key)
                if item is not None:
                    return item

    def __getitem__(self, key): 
        item = self._finditem(self.config, key)

        if item is None:
            raise ConfigParamNotFound(f'Did not find parameter {key} in the config file')

        return item
