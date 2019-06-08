import os
import sys
from glob import glob
import subprocess

known_programs = {'docking': ['autodock', 'vina', 'dock', 'glide', 'moe', 'gold'], \
     'rescoring': ['autodock', 'vina', 'dock', 'glide', 'moe', 'dsx', 'colvar']}
known_programs['scoring'] = known_programs['rescoring']

single_run_scoring_programs = ['glide', 'dock']
programs_handling_ions = ['vina', 'dock']

default_minimize_options = {'charge_method': 'gas', 'ncyc': 5000, 'maxcyc': 10000, 'cut': 999.0, 'solvent': 'vacuo'}

class ConfigSetup(object):

    def __init__(self, task, config):

        self.task = task
        self.section = task.upper()

        self.setup_instances(task, config)
        self.set_site_options(config)

    def setup_instances(self, task, config):
        self.instances = []

        if config.has_option(self.section, 'program'):

            instances = config.get(self.section, 'program').lower()
            instances = map(str.strip, instances.split(','))

            for instance in instances:
                program = ''.join([c for c in instance if not c.isdigit()]) # get program's exact name
                if program not in known_programs[task]:
                    raise ValueError("%s programs should be one of "%task.capitalize() + ", ".join(known_programs[task]))
                sys.modules[program] = __import__('dockbox.'+program, fromlist=['a'])

                options = {}
                # check if all needed executables are available
                if hasattr(sys.modules[program], 'required_programs'):
                    required_programs = getattr(sys.modules[program], 'required_programs')
                    for exe in required_programs:
                        try:
                            subprocess.check_call('which %s > /dev/null'%exe, shell=True)
                        except subprocess.CalledProcessError:
                            raise ValueError('Executable %s needed for docking with %s not found! \
Make sure the program has been installed and is in your PATH!'%(exe, program))

                # check if mandatory options are set up
                if hasattr(sys.modules[program], 'mandatory_settings'):
                    madatory_settings = getattr(sys.modules[program], 'mandatory_settings')
                    config_d = dict(config.items(instance.upper()))
                    for setting in madatory_settings:
                        if setting not in config_d or not config_d[setting]:
                            raise ValueError('Option %s when using %s is mandatory!'%(setting,program))

                # load default parameters
                if hasattr(sys.modules[program], 'default_settings'):
                    default_settings = getattr(sys.modules[program], 'default_settings')
                    for key, value in default_settings.iteritems():
                        options[key] = value

                known_settings = {}
                if hasattr(sys.modules[program], 'known_settings'):
                    known_settings = getattr(sys.modules[program], 'known_settings')

                def check_value(key, value, instance):
                   if not key in default_settings.keys():
                       raise ValueError("Option %s not recognized in instance %s!"%(key, instance))
                   # TODO: check that value has the required type, e.g. set known_settings as a dict with the type and the list of possible choices if any!
                   if key in known_settings:
                       for known_value in known_settings[key]:
                           if value.lower() == known_value.lower():
                               return known_value
                       raise ValueError("Value %s not recognized for option %s in instance %s!"%(value, key, instance))
                   elif key.endswith('dir'): # path value
                       return os.path.abspath(value)
                   else:
                       return value

                # get parameters from config file (would possibly overwrite default preset parameters)
                if config.has_section(instance.upper()):
                   config_d = dict(config.items(instance.upper()))
                   for key, value in config_d.iteritems(): 
                       options[key] = check_value(key, value, instance)

                self.instances.append((instance, program, options))
        else:
            raise ValueError("option program in section %s is required in config file!"%self.section)

    def set_site_options(self, config):
        """set options for the binding site"""

        site = {}
        required_options = ['center', 'boxsize']

        if config.has_option('DOCKING', 'site'):
            sitenames = config.get('DOCKING', 'site').lower()
            sitenames = map(str.strip, sitenames.split(','))
            for idx, name in enumerate(sitenames):
                site['site'+str(idx+1)] = [name]
                for option in required_options:
                    section = name.upper()
                    if config.has_option(section, option):
                        value = config.get(section, option)
                        site['site'+str(idx+1)].append(value)
                    else:
                        raise ValueError("Option %s in section %s is required in config file!"%(option,section))
        else:
            section = 'SITE'
            site['site1'] = [None]
            for option in required_options:
                if config.has_option(section, option):
                    value = config.get(section, option)
                    site['site1'].append(value)
                else:
                    raise ValueError("Option %s in section %s is required in config file for local docking!"%(option,section))
        self.site = site
        self.nsites = len(site)


    def get_value_yesno_option(self, config, section, option, default=False):

        if config.has_option(section, option):
            yesno = config.get(section, option).lower()
            if yesno == 'yes':
                return True
            elif yesno == 'no':
                return False
            else:
                raise ValueError("option %s should be yes or no!"%option)
        else:
            return default


    def get_value_cleanup_option(self, config, section, default=0):

        if config.has_option(section, 'cleanup'):
            value = config.get(section, 'cleanup').lower()
            if value == 'no' or value == '0':
                return 0
            elif value == 'yes' or value == '1':
                return 1
            elif value == '2':
                return 2
            else:
                raise ValueError("cleanup option in section DOCKING should be yes, no or 0, 1 or 2!")
        else:
            return default

class DockingSetup(ConfigSetup):

    def __init__(self, config):

        super(DockingSetup, self).__init__('docking', config)

        self.cleanup = self.get_value_cleanup_option(config, 'DOCKING')
        self.minimize = self.set_minimization_options(config)

        if config.has_option('DOCKING', 'clustering'):
            self.cutoff_clustering = config.getfloat('DOCKING', 'clustering')
        else:
            self.cutoff_clustering = 0.0

    def set_minimization_options(self, config):
        """set options for minimization"""

        self.minimize_options = {}
        self.minimize_options['minimization'] = self.get_value_yesno_option(config, 'DOCKING', 'minimize')

        section = 'MINIMIZATION'
        if self.minimize_options['minimization']:

            # check AMBER version
            self.minimize_options['amber_version'] = self.check_amber_version()

            # load default parameters
            for key, value in default_minimize_options.iteritems():
                self.minimize_options[key] = value

            # get parameters from config file (would possibly overwrite default preset parameters)
            if config.has_section(section):
              config_m = dict(config.items(section))
              for key, value in config_m.iteritems():
                  self.minimize_options[key] = value

        return self.minimize_options

    def check_amber_version(self):
        error_msg = 'AmberTools serial version >= 14 and <= 17 is required for minimization with DockBox!'

        if os.environ.get('AMBERHOME'):
            for exe in ['tleap', 'sander', 'cpptraj']:
                try:
                    subprocess.check_call('which %s > /dev/null'%exe, shell=True)
                except subprocess.CalledProcessError:
                    raise ValueError('Executable %s is not found in your PATH! %s'%(exe, error_msg))

            docfile = glob(os.environ.get('AMBERHOME')+'/doc/Amber*.pdf')
            amber_version = os.path.basename(docfile[0])[5:-4]
            try:
                int(amber_version)
                if amber_version not in ['14', '15', '16', '17']:
                    raise ValueError("Amber version %s detected! %s"%error_msg)
                return amber_version
            except ValueError:
                raise ValueError("Amber version not detected! %s"%error_msg)
        else:
            raise ValueError("AMBERHOME is not set! %s"%error_msg)

class RescoringSetup(ConfigSetup):

    def __init__(self, config):
        self.is_rescoring = self.get_value_yesno_option(config, 'DOCKING', 'rescoring')

        if self.is_rescoring:
            super(RescoringSetup, self).__init__('rescoring', config)

class ScoringSetup(ConfigSetup):
    pass
