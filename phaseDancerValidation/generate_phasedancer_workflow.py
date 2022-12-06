import yaml

def load_yaml(yaml_fn):
    with open(yaml_fn, 'r') as f:
        return yaml.safe_load(f)
    
def save_yaml(yaml_fn, data):
    with open(yaml_fn, 'w') as f:
        yaml.dump(data, f, default_flow_style = False)

from collections import namedtuple

sample_namedtuple_fields = ['chemistry', 'accuracy', 'coverage', 'ttype', 'region', 'sim_number', 'name', 'sample', 'simulation_dir', 'phasedancer_dir', 'phasedancer_bin', 'end']

Sample = namedtuple('Sample', sample_namedtuple_fields )
Contig = namedtuple('Contig', sample_namedtuple_fields + ['contig'])


def copying_reference(sample):
    
    text = '# copying reference\n\n'
    #text += 'mkdir -p {phasedancer_dir}data/{sample}\n'.format(**sample._asdict())
    simulated_reads = '{simulation_dir}data/simulations/{name}/sim{sim_number}/mutatedseq-{region}-{name}-sim{sim_number}-all.fasta'.format(**sample._asdict())
    text += 'cp ' + simulated_reads + ' {phasedancer_dir}data//mutatedseq-{region}-{name}-sim{sim_number}-all.fasta\n\n'.format(**sample._asdict())
    
    return text

def placing_sample(sample):
    
    text = '# placing sample\n\n'
    text += 'mkdir -p {phasedancer_dir}data/{sample}\n'.format(**sample._asdict())
    simulated_reads = '{simulation_dir}data/simulations/{name}/sim{sim_number}/{region}-{name}-sim{sim_number}-{ttype}-{chemistry}-{accuracy}-{coverage}x/{region}-{name}-sim{sim_number}-{ttype}-{chemistry}-{accuracy}-{coverage}x.fastq'.format(**sample._asdict())
    text += 'cp ' + simulated_reads + ' {phasedancer_dir}data/{sample}/{sample}.fastq\n\n'.format(**sample._asdict())
    
    return text
    
def placing_contig(contig):
    
    text = '# placing contig\n\n'
    text += 'mkdir -p {phasedancer_dir}data/{sample}/{contig}-start-mutatedseq-{region}-{name}-sim{sim_number}\n'.format(**contig._asdict())
    start_contig = '{simulation_dir}data/simulations/{name}/sim{sim_number}/starts/{contig}-start-mutatedseq-{region}-{name}-sim{sim_number}.fasta'.format(**contig._asdict())
    text += 'cp ' + start_contig + ' {phasedancer_dir}data/{sample}/{contig}-start-mutatedseq-{region}-{name}-sim{sim_number}/{contig}-start-mutatedseq-{region}-{name}-sim{sim_number}.fasta\n\n'.format(**contig._asdict())
    
    return text
    
def creating_index(sample):        
    
    text = '# creating index\n\n'
    text += '(cd {phasedancer_bin} && snakemake --snakefile=GenerateIndices.smk --config sample={sample} --cores 5)\n\n'.format(**sample._asdict())
    
    return text
    
def loading_index(sample):
    
    text = '# loading index\n\n'
    text += '{phasedancer_bin}load_index.sh {sample} 1\n\n'.format(**sample._asdict())
    
    return text
    
def starting_snakemake(phasedancer_bin):
    
    text = '# starting snakemake\n\n'
    text += f'(cd {phasedancer_bin} && snakemake --config sample=all --cores 80 --snakefile=Validate.smk --rerun-incomplete) \n\n'  
    
    return text
            
def unloading_index(sample):
    
    text = '# unloading index\n\n'
    text += '{phasedancer_bin}unload_index.sh {sample} 1\n\n'.format(**sample._asdict())        
    
    return text


def copying_contig(sample):
    
    text = '# copying contig\n\n'
    text += 'cp {phasedancer_dir}data/{region}-{name}-sim{sim_number}-{ttype}-{chemistry}-{accuracy}-{coverage}x/{contig}-start-mutatedseq-{region}-{name}-sim{sim_number}/minimap2/nc_000/merged_contigs/seq_000.seq_*.{contig}-start*'
    text += ' {simulation_dir}data/simulations/{name}/sim{sim_number}/{region}-{name}-sim{sim_number}-{ttype}-{chemistry}-{accuracy}-{coverage}x/{contig}.phaseDancer.fasta\n\n'       
    
    return text.format(**sample._asdict()) 


def generate_script(config_fn, simulation_dir, phasedancer_dir, phasedancer_bin):
    
    config = load_yaml(config_fn)
    
    samples = []
    contigs = []
    
    pd_config = { 'samples': {} }
    
    text = '# setting environment variables\n\n'
    
    print(f'export PHASEDANCER_SIMULATOR_DATA_DIR={simulation_dir}')
    print(f'export PHASEDANCER_DATA_DIR={phasedancer_dir}')
    print(f'export PHASEDANCER_BIN={phasedancer_bin}\n\n')
    
    
    for name, simulation in config['simulations'].items():
        if 'start-contig' not in simulation:
            continue
     
        ttype = simulation['type']
        contig_size = simulation['start-contig']['contig-size']
        contig_extension_size = simulation['start-contig']['contig-extension-size']
        iterations_right = simulation['start-contig']['iterations']
        
        for chemistry in simulation['chemistries']:
     
            for accuracy in simulation['accuracies']:
               
                for coverage in simulation['coverages']:
                    #if coverage != 20: continue
                    
                    for region in simulation['regions']:
                        for sim_number in range(simulation['simulations-number']):

                            sample_name = f"{region}-{name}-sim{sim_number}-{ttype}-{chemistry}-{accuracy}-{coverage}x"
                            context = {'chemistry': chemistry, 'accuracy': accuracy, 'coverage': coverage, 'ttype': ttype, 'region': region, 'sim_number': sim_number,'sample': sample_name, 'name': name, 'simulation_dir': simulation_dir, 'phasedancer_dir': phasedancer_dir, 'simulation_dir': simulation_dir, 'phasedancer_bin': phasedancer_bin, 'end':  iterations_right}

                            sample = Sample(**context)
                            samples.append(sample)

                            pd_config['samples'][sample_name] = {}
                            pd_sample_config = pd_config['samples'][sample_name] 

                            pd_sample_config['technology'] = 'pb' if chemistry in ['P4C2', 'P5C3', 'P6C4'] else 'ont'
                            pd_sample_config['contig-size'] = contig_size
                            pd_sample_config['contig-extension-size'] = contig_extension_size 
                            pd_sample_config['iterations'] = iterations_right
                            
                            contig_names_extended = []

                            for contig_name in simulation['start-contig']['contigs']:
                                contig = Contig(**context, contig = contig_name)

                                contig_name_extended = f'{contig_name}-start-mutatedseq-{region}-{name}-sim{sim_number}'
                                contig_names_extended.append(contig_name_extended)

                                contigs.append(contig)

                            pd_sample_config['contigs'] = contig_names_extended
     

    print(''.join(map(placing_sample, samples)))
    print(''.join(map(placing_contig, contigs)))
   
    print(''.join(map(creating_index, samples)))
    print(''.join(map(loading_index, samples)))
    print(starting_snakemake(phasedancer_bin))
    print(''.join(map(unloading_index, samples)))
    print(''.join(map(copying_contig, contigs)))

    
    save_yaml('phaseDancer_config.yaml', pd_config)
 

generate_script('config-flat-2.yaml',  ' /home/basia/simulator/PhaseDancerSimulator/', ' /home/basia/simulator/assemblies/', '/home/basia/phaseDancer2/')

