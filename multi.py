import sys
import subprocess

known_program_names = {'docking': ['autodock', 'vina', 'dock', 'glide', 'moe'], 'rescoring': ['vina', 'autodock', 'glide']}

class MultiProgramTask(object):

    def __init__(self, task, config):
        self.task = task
        self.section = task.upper()

        if self.task in known_program_names:
            self.known_program_names = known_program_names[self.task]
        else:
            raise ValueError("argument task should be one of " + ", ".join(known_program_names.keys()))

        self.setup_instances(task, config)
        self.set_site_options(config)

    def setup_instances(self, task, config):
        self.instances = []

        if config.has_option(self.section, 'program'):

            instances = config.get(self.section, 'program').lower()
            instances = map(str.strip, instances.split(','))

            for instance in instances:
                program = ''.join([c for c in instance if not c.isdigit()]) # get program's exact name
                if program not in self.known_program_names:
                    raise ValueError("%s programs should be one of "%task.capitalize() + ", ".join(self.known_program_names))
                __import__(program)

                options = {}
                # check if all needed executables are available
                if hasattr(sys.modules[program], 'required_programs'):
                    required_programs = getattr(sys.modules[program], 'required_programs')
                    for exe in required_programs:
                        try:
                            subprocess.check_call('which %s > /dev/null'%exe, shell=True)
                        except subprocess.CalledProcessError:
                            raise ValueError('Executable %s needed for docking with %s is not found in your PATH! \
Make sure the program has been installed!'%(exe,program))

                # load default parameters
                if hasattr(sys.modules[program], 'default_settings'):
                    default_settings = getattr(sys.modules[program], 'default_settings')
                    for key, value in default_settings.iteritems():
                        options[key] = value

                # check config file (would possibly default preset parameters)
                if config.has_section(instance.upper()):
                   config_d = dict(config.items(instance.upper()))
                   for key, value in config_d.iteritems():
                       options[key] = value

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

    def is_yesno_option(self, config, section, option, default=False):
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


class MultiProgramScoring(MultiProgramTask):

    def __init__(self, config):
        super(MultiProgramScoring, self).__init__('rescoring', config)

class MultiProgramDocking(MultiProgramTask):

    def __init__(self, config):
        super(MultiProgramDocking, self).__init__('docking', config)

        self.cleanup = self.is_yesno_option(config, 'DOCKING', 'cleanup')
        self.minimize = self.is_yesno_option(config, 'DOCKING', 'minimize')
