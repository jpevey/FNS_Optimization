# Average_Run_2
# This mess is brought to you by John Pevey
import random
from deap import base
from deap import creator
from deap import tools
import collections
import os
import math
import time

# This value sets the number of zones, or values to solve for.
# 1 - Fast and Thermal zone have same radius. No plates
# 2 - Fast and thermal have different values. No plates
# 3 - Fast and thermal have different values. Added cad and uranium plate
# 4 - Fast 1, fast 2, and thermal zones have different radius values, with both plates.
# 5 - Fast 1, fast 2, and thermal zones have different radius values, no plates
job_type = 3

# Random seed for python run
random.seed(97654356787654)
# pick job type 8inch or 10inch
job_zone_type = "10inch"
debug = False
# restarting function flag and filename
restart_from_file = False
restart_file_string = "Restart_File.txt"

#Types of possible mutation schemes:
# random_mutation - 1-all variables mutate to any value in possible range
# balanced_mass - when fuel radius mutates, other fuel mutates to keep fuel constant
# balanced_mutation - when a fuel or plate mutates, other fuel or plate mutate opposite way
mutation_type = 'balanced_mutation'
# If True, it does the classic crossover. Two parents are combiined to give an offspring, i.e. 1 from parent 1 3 from parent 2, etc...
# types: 'beta_crossover' 'classic_crossover' 'averaging_crossover'
crossover_type = 'classic_crossover'
#Maximum value for beta crossover  is 1.0 + beta_delta
beta_delta = 0.3

#Significant digits of fuel radius and plate thickness
sig_dig_fuel = 2
sig_dig_plate = 2
# Maximum value that a plate can have:
maximum_plate_value = 0.95
minimum_plate_value = 0.3
minimum_plate_thickness = 0.1
full_clear_of_files = True
make_Dir = False
Run_Directory = os.getcwd()
Hide_Mcnp = True
run_MCNP_bool = True
file_count = 0
generation = 0
run_on_necluster = True


children_list = []
use_three_zones = False

cluster_input_string = """#!/bin/bash

#PBS -V
#PBS -q fill
#PBS -l nodes=1:ppn=8

hostname
module unload mpi
module load intel/12.1.6
module load openmpi/1.6.5-intel-12.1
module load MCNP6

RTP="/tmp/runtp--".`date "+%R%N"`
cd $PBS_O_WORKDIR
mcnp6.mpi TASKS 8 name=%%%INPUT%%% runtpe=$RTP
grep -a "final result" %%%INPUT%%%o > %%%INPUT%%%_done.dat
rm $RTP"""

# Chance of mutation in new child
MUTPB = .40
population_int = 100
number_of_generations = 100
keepers_int = 20
keff_target = 0.95
mass_target = 2500.0
# The +/- change on random variables in crossover when using averaging method, i.e. use_classic_crossover = false
random_window = 0.25
# The +/- change on random variables in mutations
random_mutation_window = 0.90
# keff range to run the true source calculation
keff_min = 0.9495
keff_max = 0.9505
# minimum number of source jobs to run per generation. the keff bins expand to get this many source jobs
number_of_source_jobs_to_run = 100
maximum_keff_of_source_job = 0.9525
best_flux_val = 0.0

if job_zone_type == "10inch":
    template_file_str = "Template_File_10inch_Stage_"+str(job_type)+".txt"
    if use_three_zones == True:
        template_file_str = "fixthis"
if job_zone_type == "8inch":
    print("8inch flavor of this is broken. Related to uranium/cad/void in template and in buildmcnp...")
    exit(0)
model_dict = collections.OrderedDict()
run_command = "mcnp6 tasks 8 inp="
keff_run_str = "kcode 5000 1.000000 20 120 15000\nksrc  36.261700 23.528000 39.365000"
source_run_str = "sdef pos 51.6365 22.55 37.6747 erg=2.5 par=1\nNPS 1e5"
if job_zone_type == "10inch":
    source_run_str = "sdef pos 61.7965 22.55 37.6747 erg=2.5 par=1\nNPS 1e5"

# log file
output_file_str = "Output.csv"

# model variables
cadmium_density = -8.65  # g/cc (negative for mcnp)
uranium_density = -19.1  # g/cc (negative for mcnp)
clad_thickness_cm = 0.1  # cm
maximum_fuel_radius = 2.54 / 2 - clad_thickness_cm - .001
minimum_fuel_radius = 0.1
print(maximum_fuel_radius)

pin_length = (8 - (clad_thickness_cm / 2.54) * 2) * 2.54  # in cm
if job_zone_type == "10inch":
    pin_length = (10 - (clad_thickness_cm / 2.54) * 2) * 2.54  # in cm
IND_SIZE = 1


# clears
def full_Mcnp_Run_Cleanup(Run_bool, pass_list):
    if Run_bool == True:
        list_of_files = os.listdir()

        # print(list_of_files)
        for file in list_of_files:
            pass_on_file = False
            if '_ran_as_input' in file:
                continue

            for pass_file in pass_list:
                if pass_file in file:
                    print("Passing on %s, since it's a parent of next generation" % (file))
                    pass_on_file = True
            if pass_on_file == True:
                continue
            if '.inp' in file:
                os.remove(file)
            if '.out' in file:
                try:
                    os.remove(file)
                except:
                    pass
            if 'runt' in file:
                os.remove(file)
            if 'src' in file:
                os.remove(file)
            if 'mesh' in file:
                os.remove(file)


full_Mcnp_Run_Cleanup(full_clear_of_files, [])


def light_mcnp_run_cleanup(Run_bool, pass_list):
    if Run_bool == True:
        list_of_files = os.listdir()

        # print(list_of_files)
        for file in list_of_files:
            pass_on_file = False
            if '_ran_as_input' in file:
                continue

            for pass_file in pass_list:
                if pass_file in file:
                    print("Passing on %s, since it's a parent of next generation" % (file))
                    pass_on_file = True
            if pass_on_file == True:
                continue
            if 'runt' in file:
                os.remove(file)
            if 'src' in file:
                os.remove(file)


# Setting up variables for each MCNP run. i.e, setting what will be maximized and what will be minimized?
# weights = k effective, Total_Mass,  Flux in Experiment Zone
creator.create("FitnessMulti", base.Fitness, weights=(1.0, 1.0, 1.0, 1.0))
creator.create("Individual", list, fitness=creator.FitnessMulti)

toolbox = base.Toolbox()


class Assembly:
    def __init__(self):
        # Initial Model Variables
        self.fast_fuel_radius = round(random.uniform(minimum_fuel_radius, maximum_fuel_radius), sig_dig_fuel)
        self.fast_fuel_2_radius = round(random.uniform(minimum_fuel_radius, maximum_fuel_radius), sig_dig_fuel)
        self.thermal_fuel_radius = round(random.uniform(minimum_fuel_radius, maximum_fuel_radius), sig_dig_fuel)
        # This value is the split between the void volume and the cadmium/uranium plater
        self.void_value = round(random.uniform(minimum_plate_value, maximum_plate_value), sig_dig_plate)
        # This value is the split between the cadmium and the uranium plates (given the void split above)
        self.cad_u_ratio = round(random.uniform(minimum_plate_value, maximum_plate_value), sig_dig_plate)

        # Job with just 1 fuel radius
        if job_type == 1:
            self.thermal_fuel_radius = 0
            self.cad_u_ratio = 1.0
            self.void_value = 1.0
            self.fast_fuel_2_radius = 0
        # Job with fast and thermal fuel radius
        if job_type == 2:
            self.cad_u_ratio = 1.0
            self.void_value = 1.0
            self.fast_fuel_2_radius = 0
        # Job with fast and thermal fuel radius and plates
        if job_type == 3:
            self.fast_fuel_2_radius = 0

        self.input_name = "Null"
        self.run_keff = ""
        self.keff = "Null"
        self.fast_flux = 0.0
        self.thermal_flux = 0.0

        self.keff_score = 'Null'
        self.generation = 1
        self.mass = 'Null'
        self.mass_score = 'Null'
        self.ran_true_source = 'Null'
        self.is_parent = False


def job_name():
    global file_count
    global generation

    file_count = file_count + 1

    return ("input_" + str(file_count) + "_gen_" + str(generation))


toolbox.register("Assembly", Assembly)
toolbox.register("job_name", job_name)
toolbox.register("individual", tools.initCycle, creator.Individual,
                 (toolbox.Assembly,
                  toolbox.job_name
                  ),
                 n=IND_SIZE)
toolbox.register("population", tools.initRepeat, list, toolbox.individual)


def get_keff(output_file_string):
    output_file = open(output_file_string, 'r')
    found_keff = False
    for line in output_file:
        if "the final estimated combined collision/absorption/track-length keff =" in line:
            line_split_1 = line.split('the final estimated combined collision/absorption/track-length keff = ')
            line_split_2 = line_split_1[1].split(' with ')
            keff = line_split_2[0]
            print(keff)
            found_keff = True
        if " surface  2    " in line:
            get_tally = True
            continue
    if found_keff == False:
        print("Unable to find a keff for input " + output_file_string)
        return 0.0
    return keff


def score_keff(keff):
    # Calculating the score. If the keff is closer than +/-.0025, it gets a max score
    try:
        keff_score = 1 / abs((float(keff_target) - float(keff)))
    except(ZeroDivisionError):
        keff_score = 4.0
    if keff_score > 400:
        keff_score = 400
    keff_score = keff_score / 100
    print(keff_score)

    return keff_score, keff


def check_duplicate(pop):
    print("Checking for duplicates in this generation...")
    count = 0
    ind_count = 0
    for ind in pop:
        ind_count += 1
        # skips the parents
        if count < keepers_int:
            count += 1
            continue
        ind_2_count = 0
        for ind_2 in pop:
            ind_2_count += 1
            # skips itself...
            if ind_2_count == ind_count:
                continue
            if ind[0].fast_fuel_radius == ind_2[0].fast_fuel_radius:
                if ind[0].thermal_fuel_radius == ind_2[0].thermal_fuel_radius:
                    if ind[0].void_value == ind_2[0].void_value:
                        if ind[0].cad_u_ratio == ind_2[0].cad_u_ratio:
                            # Job with just 1 fuel radius
                            if job_type == 1:
                                random_int = 1
                            # Job with fast and thermal fuel radius
                            if job_type == 2:
                                random_int = random.randint(1, 2)
                            # Job with fast and thermal fuel radius and plates
                            if job_type == 3:
                                random_int = random.randint(1, 4)
                            # Job with fast1, fast2, thermal, and both plates
                            if job_type == 4:
                                random_int = random.randint(1, 5)
                            print("Found a duplicate! #", ind_count, "and #", ind_2_count,
                                  "are the same. Changing val:", random_int)
                            if random_int == 1:
                                ind[0].fast_fuel_radius = get_new_value(ind[0].fast_fuel_radius,'fast_fuel_radius')
                            if random_int == 2:
                                ind[0].thermal_fuel_radius = get_new_value(ind[0].thermal_fuel_radius,'thermal_fuel_radius')
                            if random_int == 3:
                                ind[0].cad_u_ratio = get_new_value(ind[0].cad_u_ratio,'cad_u_ratio')
                            if random_int == 4:
                                ind[0].void_value = get_new_value(ind[0].void_value,'void_value')
                            if random_int == 5:
                                ind[0].fast_fuel_2_radius = get_new_value(ind[0].fast_fuel_2_radius,'fast_fuel_2_radius')
        count += 1
    return pop


def calculate_mass(individual):
    fast_fuel_radius = float(individual[0].fast_fuel_radius)
    thermal_fuel_radius = float(individual[0].thermal_fuel_radius)
    cadmiun_uranium_mix = float(individual[0].cad_u_ratio)

    # fast zone mass calc
    number_of_fast_fuel_pins = 8 * 18 + 9 + 9 * 18
    fz_volume = math.pi * fast_fuel_radius ** 2 * pin_length * number_of_fast_fuel_pins
    fz_mass_kg = fz_volume * (-1 * uranium_density) / 1000

    if use_three_zones == True:
        # calculating for zone 1 pins
        number_of_fast_fuel_pins = 8 * 18 + 9
        fz_volume = math.pi * fast_fuel_radius ** 2 * pin_length * number_of_fast_fuel_pins

        # calculating for zone 2 pins
        number_of_fast_fuel_pins = 9 * 18
        # Combined volume of zone 1 and 2
        fz_volume = fz_volume + math.pi * fast_fuel_radius ** 2 * pin_length * number_of_fast_fuel_pins
        fz_mass_kg = fz_volume * (-1 * uranium_density) / 1000

    # thermal xone mass calc
    number_of_thermal_fuel_pins = 8 * 18
    tz_volume = math.pi * thermal_fuel_radius ** 2 * pin_length * number_of_thermal_fuel_pins
    tz_mass_kg = tz_volume * (-1 * uranium_density) / 1000

    # fission plate mass calc
    # mass of plate in kg

    void_depth = individual[0].void_value * 1.27
    cadmium_depth = (1.27 - void_depth) * individual[0].cad_u_ratio
    uranium_depth = (1.27 - void_depth) * (1 - individual[0].cad_u_ratio)

    uranium_plate_volume_cc = uranium_depth * 87.5030 * 95.8850
    uranium_plate_mass_g = (-1) * uranium_density * uranium_plate_volume_cc
    uranium_plate_mass_kg = uranium_plate_mass_g / 1000
    print(uranium_depth,uranium_plate_volume_cc,uranium_plate_mass_g,uranium_plate_mass_kg)
    total_mass_kg = fz_mass_kg + tz_mass_kg + uranium_plate_mass_kg

    total_mass_lb = total_mass_kg * 2.20462

    return total_mass_lb

def get_fast_and_thermal(output_file_string):
    in_tally = False
    flux_value = 0.0
    thermal_flux_value = 0.0
    with open(output_file_string, 'r') as output_file:
        for line in output_file:
            line = line.strip()
            if line.startswith("1tally        4"):
                in_tally = True
            if in_tally == True:
                if line.startswith("1.0000E+01   "):
                    line_split = line.split(' ')
                    flux_value = line_split[3]
                if line.startswith("1.0000E-02   "):
                    line_split = line.split(' ')
                    thermal_flux_value = line_split[3]
                if line.startswith("total   "):
                    in_tally = False
    return flux_value, thermal_flux_value

# This his how fast flux is scored. Previously, all scores above 0.001 were given max score. Now, it's a tiered approach
# with each +0.0001 giving a higher score
def score_fast_flux(flux_value):
    flux_value = float(flux_value)
    flux_score = round(flux_value, 4)

    return flux_score

def score_mass(mass):
    if mass <= mass_target:
        return 1.0
    if mass > mass_target:
        if (1 - (mass - mass_target) * 0.001) < 0:
            return 0.
        return (1 - (mass - mass_target) * 0.001)


def evaluate(individual):
    fast_flux_score = 0.
    thermal_flux_score = 0.
    keff_score = 0.

    # getting variables, naming specifically

    fast_fuel_radius = individual[0].fast_fuel_radius
    thermal_fuel_radius = individual[0].thermal_fuel_radius
    cadmiun_uranium_mix = individual[0].cad_u_ratio
    input_filename = individual[0].input_name

    if run_on_necluster == True:
        keff_score, keff = score_keff(get_keff(input_filename + "_keff.inpo"))
    if run_on_necluster == False:
        keff_score, keff = score_keff(get_keff(input_filename + "_keff.inp.out"))
    individual[0].keff = keff
    individual[0].mass = calculate_mass(individual)
    individual[0].keff_score = keff_score
    mass_score = score_mass(individual[0].mass)
    individual[0].mass_score = mass_score
    if debug:
        print("Eval file: %s w/ vars: %r, %r, %r. Keff Score: %s Mass Score: %s" % (
            input_filename, fast_fuel_radius, thermal_fuel_radius, cadmiun_uranium_mix, keff_score, mass_score))

    # will evaluate the total fast flux i
    if individual[0].ran_true_source == True:
        if run_on_necluster == True:
            fast_flux_value, thermal_flux_value = get_fast_and_thermal(input_filename + "_source.inpo")
        if run_on_necluster == False:
            fast_flux_value, thermal_flux_value = get_fast_and_thermal(input_filename + "_source.inp.out")
        fast_flux_score = score_fast_flux(fast_flux_value)
        thermal_flux_score = 1 - float(thermal_flux_value)
        individual[0].fast_flux = fast_flux_value
        individual[0].thermal_flux = thermal_flux_value

    print(keff_score, mass_score, fast_flux_score, thermal_flux_score)
    return float(keff_score), float(mass_score), float(fast_flux_score), float(thermal_flux_score)


# registering the above function as evaluate
toolbox.register("evaluate", evaluate)


# creates mcnp input from template and model params
def build_mcnp_model(individual):
    print("Writing out MCNP model")

    if individual[0].is_parent == True:
        return individual

    fast_fuel_radius = str(individual[0].fast_fuel_radius)
    if fast_fuel_radius == maximum_fuel_radius:
        fast_fuel_radius = maximum_fuel_radius - 0.0000001
    fast_fuel_clad_radius = str(individual[0].fast_fuel_radius + clad_thickness_cm)
    fast_fuel_2_radius = str(individual[0].fast_fuel_2_radius)
    if fast_fuel_2_radius == maximum_fuel_radius:
        fast_fuel_2_radius = maximum_fuel_radius - 0.0000001
    fast_fuel_2_clad_radius = str(individual[0].fast_fuel_2_radius + clad_thickness_cm)
    thermal_fuel_radius = str(individual[0].thermal_fuel_radius)
    if thermal_fuel_radius == maximum_fuel_radius:
        thermal_fuel_radius = maximum_fuel_radius - 0.0000001
    thermal_fuel_clad_radius = str(individual[0].thermal_fuel_radius + clad_thickness_cm)
    cadmiun_uranium_mix = str(individual[0].cad_u_ratio)

    fission_plate_split = 41.841 + (1.27) * individual[0].cad_u_ratio

    if job_zone_type == "10inch":
        void_vol_split = 52.001 + (1.27) * individual[0].void_value
        fission_plate_split = void_vol_split + (53.271 - void_vol_split) * individual[0].cad_u_ratio
        fission_plate_split = round(fission_plate_split, 5)
        if fission_plate_split >= 53.271:
            fission_plate_split = 53.270
    if void_vol_split == fission_plate_split:
        void_vol_split = void_vol_split - minimum_plate_thickness
    if individual[0].run_keff == True:
        print("Job %s has already been run in previous gen. Updating output file for next gen")

        return

    individual[0].input_name = individual[1]
    # opening template file
    template_file = open(template_file_str, 'r')

    # opening new input file for keff job
    input_keff = open(individual[0].input_name + "_keff.inp", 'w')

    for line in template_file:
        if "%fast_fuel_radius" in line:
            line = line.replace("%fast_fuel_radius", fast_fuel_radius)
        if "%fast_fuel_clad_radius" in line:
            line = line.replace("%fast_fuel_clad_radius", fast_fuel_clad_radius)
        if "%fast_fuel_2_radius" in line:
            line = line.replace("%fast_fuel_2_radius", fast_fuel_2_radius)
        if "%fast_fuel_2_clad_radius" in line:
            line = line.replace("%fast_fuel_2_clad_radius", fast_fuel_2_clad_radius)
        if "%thermal_fuel_radius" in line:
            line = line.replace("%thermal_fuel_radius", thermal_fuel_radius)
        if "%thermal_fuel_clad_radius" in line:
            line = line.replace("%thermal_fuel_clad_radius", thermal_fuel_clad_radius)
        if "%cadmiun_uranium_mix" in line:
            line = line.replace("%cadmiun_uranium_mix", fuel_coolant_material_def(cadmiun_uranium_mix, "cad"))
        if "%fission_plate_density" in line:
            line = line.replace("%fission_plate_density", uranium_cadmium_density(cadmiun_uranium_mix))
        if "%run_info" in line:
            line = line.replace("%run_info", keff_run_str)
        if "%fission_plate_split" in line:
            line = line.replace("%fission_plate_split", str(round(fission_plate_split, 5)))
        if "%void_split" in line:
            line = line.replace("%void_split", str(round(void_vol_split, 5)))
        input_keff.write(line)

    # Making input for source run
    input_source = open(individual[0].input_name + "_source.inp", 'w')
    template_file = open(template_file_str, 'r')
    for line in template_file:
        if "%fast_fuel_radius" in line:
            line = line.replace("%fast_fuel_radius", fast_fuel_radius)
        if "%fast_fuel_clad_radius" in line:
            line = line.replace("%fast_fuel_clad_radius", fast_fuel_clad_radius)
        if "%fast_fuel_2_radius" in line:
            line = line.replace("%fast_fuel_2_radius", fast_fuel_2_radius)
        if "%fast_fuel_2_clad_radius" in line:
            line = line.replace("%fast_fuel_2_clad_radius", fast_fuel_2_clad_radius)
        if "%thermal_fuel_radius" in line:
            line = line.replace("%thermal_fuel_radius", thermal_fuel_radius)
        if "%thermal_fuel_clad_radius" in line:
            line = line.replace("%thermal_fuel_clad_radius", thermal_fuel_clad_radius)
        if "%cadmiun_uranium_mix" in line:
            line = line.replace("%cadmiun_uranium_mix", fuel_coolant_material_def(cadmiun_uranium_mix, "cad"))
        if "%fission_plate_density" in line:
            line = line.replace("%fission_plate_density", uranium_cadmium_density(cadmiun_uranium_mix))
        if "%run_info" in line:
            line = line.replace("%run_info", source_run_str)
        if "%fission_plate_split" in line:
            line = line.replace("%fission_plate_split", str(round(fission_plate_split, 5)))
        if "%void_split" in line:
            line = line.replace("%void_split", str(round(void_vol_split, 5)))
        input_source.write(line)

    return individual


# given the cad/u split, returns the density
def uranium_cadmium_density(split_value):
    value = uranium_density * float(split_value) + cadmium_density * (1 - float(split_value))
    return str(value)


# takes a split value and returns a string with a mcnp material definitions
def fuel_coolant_material_def(split_value, type):
    split_value = float(split_value)
    if type == 'fast':
        value_string = " 92235. " + str(-0.0975 * split_value) + "\n" + \
                       "     92238. " + str(-0.9025 * split_value) + "\n" + \
                       "     82204. " + str(-0.014 * (1 - split_value)) + "\n" + \
                       "     82206. " + str(-0.241 * (1 - split_value)) + "\n" + \
                       "     82207. " + str(-0.221 * (1 - split_value)) + "\n" + \
                       "     82208. " + str(-0.524 * (1 - split_value))
    if type == 'thermal':
        value_string = " 92235. " + str(-0.0975 * split_value) + "\n" + \
                       "     92238. " + str(-0.9025 * split_value) + "\n" + \
                       "     1001. " + str(-0.143716 * (1 - split_value)) + "\n" + \
                       "     6000. " + str(-0.856284 * (1 - split_value))
        if 1 - split_value < 0:
            split_value = 1.0
            value_string = " 92235. " + str(-0.0975 * split_value) + "\n" + \
                           "     92238. " + str(-0.9025 * split_value)

    if type == 'cad':
        value_string = " 92235. " + str(-0.0975 * split_value) + "\n" + \
                       "     92238. " + str(-0.9025 * split_value) + "\n" + \
                       "     48000. " + str(-1.0 * (1 - split_value))
        if 1 - split_value < 0:
            split_value = 1.0
            value_string = " 92235. " + str(-0.0975 * split_value) + "\n" + \
                           "     92238. " + str(-0.9025 * split_value)
    return value_string


# runs the current generation of jobs
def run_Mcnp(variable, generation, run_type):
    mcnp_input_file = "input_"

    generation_string = "_gen_" + str(generation)
    ##print(os.getcwd())
    os.chdir(Run_Directory)
    ##print(os.getcwd())

    list_of_files = os.listdir()

    # print(list_of_files)

    for file in list_of_files:
        if mcnp_input_file in file:

            if '_ran_as_input' in file:
                print("The job: " + file + " has already been run, continuing")
                continue

            if '.out' in file:
                continue
            if variable not in file:
                continue
            if run_type not in file:
                continue
            if generation_string not in file:
                continue

            # print("Running MCNP for file:" + file)
            print("Running/Submitting MCNP job: " + file)
            if Hide_Mcnp == True:
                if run_MCNP_bool:
                    os.system(run_command + file + " out=" + file + ".out > file_name.whatever")
                    light_mcnp_run_cleanup(full_clear_of_files, [])
                pass
            if Hide_Mcnp == False:
                if run_MCNP_bool:
                    os.system(run_command + file + " out=" + file)
                    light_mcnp_run_cleanup(full_clear_of_files, [])
                pass
            os.chdir(Run_Directory)


def run_on_cluster(variable, generation, run_type):
    mcnp_input_file = "input_"

    generation_string = "_gen_" + str(generation)
    ##print(os.getcwd())
    os.chdir(Run_Directory)
    ##print(os.getcwd())

    list_of_files = os.listdir()

    list_of_jobs_submitted = []
    # print(list_of_files)

    for file in list_of_files:
        if mcnp_input_file in file:
            if '_ran_as_input' in file:
                print("The job: " + file + " has already been run, continuing")
                continue
            if '.out' in file:
                continue
            if variable not in file:
                continue
            if run_type not in file:
                continue
            if generation_string not in file:
                continue
            print("Creating cluster script for  MCNP job: " + file)

            script_string = cluster_input_string
            script_string = script_string.replace("%%%INPUT%%%", file)
            script_file_string = file + "_script.txt"
            script_file = open(script_file_string, 'w')
            script_file.write(script_string)
            script_file.close()

            print("Submitting MCNP job: " + file)
            list_of_jobs_submitted.append(file)

            os.system('qsub ' + script_file_string)
    return list_of_jobs_submitted


toolbox.register("build_mcnp_model", build_mcnp_model)
toolbox.register("run_Mcnp", run_Mcnp)
toolbox.register("run_on_cluster", run_on_cluster)

# This function returns a new value for the given variable that's within the acceptable bounds.
def get_new_value(old_value, variable_type):
    if              variable_type == 'fast_fuel_radius' or\
                    variable_type == 'thermal_fuel_radius' or\
                    variable_type == 'fast_fuel_2_radius':
        return_variable = round(random.uniform(minimum_fuel_radius, maximum_fuel_radius), sig_dig_fuel)
        while return_variable > round((maximum_fuel_radius), 3) or \
              return_variable < round((minimum_fuel_radius), 3) :
            return_variable = maximum_fuel_radius
    if              variable_type == 'cad_u_ratio' or \
                    variable_type == 'void_value':
        return_variable = round(random.uniform(minimum_plate_value, maximum_plate_value), sig_dig_plate)

    return return_variable

def mutate_with_constant_mass(individual, variable_to_mutate):
    maximum_mass = 'null'
    minimum_mass = 'null'
    old_mass = round(calculate_mass(individual))
    new_mass = 0.0
    old_fast_val = individual[0].fast_fuel_radius
    old_therm_val = individual[0].thermal_fuel_radius

    if job_type > 3:
        include_2nd_fast_zone = True

    if variable_to_mutate == 1:
        individual[0].fast_fuel_radius = get_new_value(individual[0].fast_fuel_radius, 'fast_fuel_radius')
        new_mass = round(calculate_mass(individual))
        # Calculate the maximum change possible
        if old_mass > new_mass:
            #Calculating maximum mass we could have with thermal fuel
            individual[0].thermal_fuel_radius = maximum_fuel_radius
            maximum_mass = round(calculate_mass(individual))
            # If the maximum mass is less than the old mass, then all we can do is maximize the thermal fuel:
            if maximum_mass < old_mass:
                individual[0].thermal_fuel_radius = maximum_fuel_radius
                new_mass = round(calculate_mass(individual))
            # If the maximum mass is greater than the old_mass, then there's a value for therm radius that makes
            if maximum_mass > old_mass:
                individual[0].thermal_fuel_radius = old_therm_val
                while new_mass != old_mass:
                    individual[0].thermal_fuel_radius = individual[0].thermal_fuel_radius + 0.001
                    new_mass = round(calculate_mass(individual))
                    if new_mass > old_mass:
                        individual[0].thermal_fuel_radius = individual[0].thermal_fuel_radius - 0.001
                        new_mass = old_mass
        if old_mass < new_mass:
            individual[0].thermal_fuel_radius = minimum_fuel_radius
            minimum_mass = round(calculate_mass(individual))
            if minimum_mass > old_mass:
                individual[0].thermal_fuel_radius = minimum_fuel_radius
                new_mass = round(calculate_mass(individual))
            if minimum_mass < old_mass:
                individual[0].thermal_fuel_radius = old_therm_val
                while new_mass != old_mass:
                    individual[0].thermal_fuel_radius = individual[0].thermal_fuel_radius - 0.001
                    new_mass = round(calculate_mass(individual))
                    if new_mass < old_mass:
                        individual[0].thermal_fuel_radius = individual[0].thermal_fuel_radius + 0.001
                        new_mass = old_mass


    if variable_to_mutate == 2:
        individual[0].thermal_fuel_radius = get_new_value(individual[0].thermal_fuel_radius, 'thermal_fuel_radius')
        new_mass = round(calculate_mass(individual))

        # Calculate the maximum change possible
        if old_mass > new_mass:
            #Calculating maximum mass we could have with thermal fuel
            individual[0].fast_fuel_radius = maximum_fuel_radius
            maximum_mass = round(calculate_mass(individual))
            # If the maximum mass is less than the old mass, then all we can do is maximize the thermal fuel:
            if maximum_mass < old_mass:
                individual[0].fast_fuel_radius = maximum_fuel_radius
                new_mass = round(calculate_mass(individual))
            # If the maximum mass is greater than the old_mass, then there's a value for therm radius that makes
            if maximum_mass > old_mass:
                individual[0].fast_fuel_radius = old_fast_val
                while new_mass != old_mass:
                    individual[0].fast_fuel_radius = individual[0].fast_fuel_radius + 0.001
                    new_mass = round(calculate_mass(individual))
                    if new_mass > old_mass:
                        individual[0].fast_fuel_radius = individual[0].fast_fuel_radius - 0.001
                        new_mass = old_mass
        if old_mass < new_mass:
            individual[0].fast_fuel_radius = minimum_fuel_radius
            minimum_mass = round(calculate_mass(individual))
            if minimum_mass > old_mass:
                individual[0].fast_fuel_radius = minimum_fuel_radius
                new_mass = round(calculate_mass(individual))
            if minimum_mass < old_mass:
                individual[0].fast_fuel_radius = old_fast_val
                while new_mass != old_mass:
                    individual[0].fast_fuel_radius = individual[0].fast_fuel_radius - 0.001
                    new_mass = round(calculate_mass(individual))
                    if new_mass < old_mass:
                        individual[0].fast_fuel_radius = individual[0].fast_fuel_radius + 0.001
                        new_mass = old_mass

    print(variable_to_mutate,"old mass, new mass", old_mass, ",", new_mass, maximum_mass, minimum_mass )



    return individual

def get_plate_thicknesses(individual):
    void_thickness = individual[0].void_value * 1.27
    cadmium_thickness = (1.27 - void_thickness) * individual[0].cad_u_ratio
    uranium_thickness = 1.27 - cadmium_thickness - void_thickness

    return void_thickness, cadmium_thickness, uranium_thickness



def get_new_plate_ratios(cad_thickness, uranium_thickness):
    print(cad_thickness,uranium_thickness)
    void_ratio = round((1.27 - cad_thickness - uranium_thickness) / 1.27, sig_dig_plate)
    cad_u_ratio = round(cad_thickness / (cad_thickness + uranium_thickness), sig_dig_plate)
    print("new void, cad ratios: ", void_ratio, cad_u_ratio)
    return void_ratio, cad_u_ratio

# mutating
def mutation(individual):
    individual[0].fast_fuel_radius = float(individual[0].fast_fuel_radius)
    individual[0].thermal_fuel_radius = float(individual[0].thermal_fuel_radius)
    individual[0].cad_u_ratio = float(individual[0].cad_u_ratio)
    individual[0].void_value = float(individual[0].void_value)
    old_0 = individual[0].fast_fuel_radius
    old_1 = individual[0].thermal_fuel_radius
    old_2 = individual[0].cad_u_ratio
    old_3 = individual[0].void_value
    if use_three_zones == True:
        old_4 = individual[0].fast_fuel_2_radius
    if mutation_type == 'random_mutation':


        individual[0].fast_fuel_radius =    get_new_value(individual[0].fast_fuel_radius, 'fast_fuel_radius')
        individual[0].thermal_fuel_radius = get_new_value(individual[0].thermal_fuel_radius, 'thermal_fuel_radius')
        individual[0].cad_u_ratio =         get_new_value(individual[0].cad_u_ratio, 'cad_u_ratio')
        individual[0].void_value =          get_new_value(individual[0].void_value, 'void_value')
        if use_three_zones == True:
            individual[0].fast_fuel_2_radius = get_new_value(individual[0].fast_fuel_2_radius, 'fast_fuel_2_radius')
        print("Mutating results: Old Vals: %s %s %s %s New Vals: %s %s %s %s" % (
            old_0, old_1, old_2, old_3, individual[0].fast_fuel_radius, individual[0].thermal_fuel_radius,
            individual[0].cad_u_ratio, individual[0].void_value))

    if mutation_type == 'balanced_mass':
        if job_type == 2:
            mutate_with_constant_mass(individual, random.randint(1, 2))
        if job_type == 3:
            mutate_with_constant_mass(individual, random.randint(1, 2))

    if mutation_type == 'balanced_mutation':
        random_val_list = []
        print("balanced mutation....")
        original_fast_fuel_val  = individual[0].fast_fuel_radius
        original_therm_fuel_val = individual[0].thermal_fuel_radius
        original_void_value = individual[0].void_value
        original_cad_u_ratio = individual[0].cad_u_ratio
        #Calculates void, cadmium, uranium thicknesses
        void_thickness, cadmium_thickness, uranium_thickness = get_plate_thicknesses(individual)
        print("original plate thicknesses v, ca, ur", void_thickness, cadmium_thickness, uranium_thickness)
        print("original ratios void, cad u", original_void_value, original_cad_u_ratio)
        if job_type == 2:
            rand_val  = random.randint(1,2)
            random_val_list = rand_val
        if job_type == 3:
            rand_val = random.randint(1, 4)
            print("Mutating", rand_val, "variables")
            for i in range(rand_val):
                value_to_change = random.randint(1, 4)

                if value_to_change not in random_val_list:
                    random_val_list.append(value_to_change)
        for rand_val in random_val_list:
            if rand_val == 1:
                print("Mutatiing fast fuel 1")
                individual[0].fast_fuel_radius = get_new_value(individual[0].fast_fuel_radius, 'fast_fuel_radius')
                if (individual[0].fast_fuel_radius - original_fast_fuel_val) > 0:
                    individual[0].thermal_fuel_radius = get_new_value(individual[0].thermal_fuel_radius,
                                                                      'thermal_fuel_radius')
                    while (individual[0].thermal_fuel_radius - original_therm_fuel_val) > 0:
                        individual[0].thermal_fuel_radius = get_new_value(original_therm_fuel_val,
                                                                          'thermal_fuel_radius')
                if (individual[0].fast_fuel_radius - original_fast_fuel_val) < 0:
                    individual[0].thermal_fuel_radius = get_new_value(original_therm_fuel_val,
                                                                      'thermal_fuel_radius')
                    while (individual[0].thermal_fuel_radius - original_therm_fuel_val) < 0:
                        individual[0].thermal_fuel_radius = get_new_value(original_therm_fuel_val,
                                                                          'thermal_fuel_radius')
            if rand_val == 2:
                print("Mutatiing thermal fuel")
                individual[0].thermal_fuel_radius = get_new_value(individual[0].thermal_fuel_radius, 'thermal_fuel_radius')
                if (individual[0].thermal_fuel_radius - original_therm_fuel_val) > 0:
                    individual[0].fast_fuel_radius = get_new_value(individual[0].fast_fuel_radius,
                                                                      'thermal_fuel_radius')
                    while (individual[0].fast_fuel_radius - original_fast_fuel_val) > 0:
                        individual[0].fast_fuel_radius = get_new_value(original_fast_fuel_val,
                                                                          'thermal_fuel_radius')
                if (individual[0].thermal_fuel_radius - original_therm_fuel_val) < 0:
                    individual[0].fast_fuel_radius = get_new_value(individual[0].fast_fuel_radius,
                                                                      'fast_fuel_radius')
                    while (individual[0].fast_fuel_radius - original_fast_fuel_val) < 0:
                        individual[0].fast_fuel_radius = get_new_value(original_fast_fuel_val,
                                                                          'fast_fuel_radius')
            # Mutating Cadmium Thickness
            if rand_val == 3:
                print("Mutating cadmium plate")
                # The new cadmium thickness is a random value between
                # the minimum plate thickness + uranium plate (unchanged) and the sum of cadmium and void,
                # minus 2 minimums to account for minimum void and uranium plates
                minimum_val = minimum_plate_value
                maximum_val = void_thickness + cadmium_thickness - 2 * minimum_plate_value
                new_cad_thickness = round(random.uniform(minimum_plate_thickness, maximum_val), sig_dig_plate)
                individual[0].void_value, individual[0].cad_u_ratio = get_new_plate_ratios(new_cad_thickness, uranium_thickness)

            # Mutating Uranium Plate Thickness
            if rand_val == 4:
                print("Mutating uranium plate")
                original_uranium_thickness = uranium_thickness
                # The new uranium thickness is a random value between
                # the minimum plate thickness + uranium plate (unchanged) and the sum of cadmium and void,
                # minus 2 minimums to account for minimum void and uranium plates
                minimum_val = minimum_plate_value
                maximum_val = void_thickness + cadmium_thickness - 2 * minimum_plate_value
                new_uranium_thickness = round(random.uniform(minimum_plate_thickness, maximum_val), sig_dig_plate)
                individual[0].void_value, individual[0].cad_u_ratio = get_new_plate_ratios(cadmium_thickness, new_uranium_thickness)

                # If the uranium plate gets thicker, the thermal fuel radius gets smaller

                if (original_uranium_thickness - new_uranium_thickness) > 0:
                    individual[0].thermal_fuel_radius = get_new_value(individual[0].thermal_fuel_radius,
                                                                      'thermal_fuel_radius')
                    while (individual[0].thermal_fuel_radius - original_therm_fuel_val) > 0:
                        individual[0].thermal_fuel_radius = get_new_value(original_therm_fuel_val,
                                                                          'thermal_fuel_radius')
                if (original_uranium_thickness - new_uranium_thickness) < 0:
                    individual[0].thermal_fuel_radius = get_new_value(original_therm_fuel_val,
                                                                      'thermal_fuel_radius')
                    while (individual[0].thermal_fuel_radius - original_therm_fuel_val) < 0:
                        individual[0].thermal_fuel_radius = get_new_value(original_therm_fuel_val,
                                                                          'thermal_fuel_radius')
        print("New fast val, old val, new therm val, old val", individual[0].fast_fuel_radius,
                  original_fast_fuel_val, individual[0].thermal_fuel_radius, original_therm_fuel_val)
    # Job with just 1 fuel radius
    if job_type == 1:
        individual[0].thermal_fuel_radius = 0
        individual[0].cad_u_ratio = 1.0
        individual[0].void_value = 1.0
        individual[0].fast_fuel_2_radius = 0
    # Job with fast and thermal fuel radius
    if job_type == 2:
        individual[0].cad_u_ratio = 1.0
        individual[0].void_value = 1.0
        individual[0].fast_fuel_2_radius = 0
    # Job with fast and thermal fuel radius and plates
    if job_type == 3:
        individual[0].fast_fuel_2_radius = 0

    return individual


pop = toolbox.population(population_int)

for ind in pop:
    #ind[0].void_value = 0.5
    #ind[0].cad_u_ratio = 0.25
    print(get_plate_thicknesses(ind))
    mutation(ind)
    print(get_plate_thicknesses(ind))



toolbox.register("mutate", mutation)

# this function selects the best
toolbox.register("select", tools.selBest)


def crossover(individual_1, individual_2, children_list):
    child_obj = Assembly()
    if crossover_type == 'beta_crossover':


        beta = random.uniform(0, 1.0 + beta_delta)
        child_obj.fast_fuel_radius    = round(beta * individual_1[0].fast_fuel_radius + (1 - beta) * individual_2[0].fast_fuel_radius, sig_dig_fuel)

        beta = random.uniform(0, 1.0 + beta_delta)
        child_obj.thermal_fuel_radius =  round(beta * individual_1[0].thermal_fuel_radius + (1 - beta) * individual_2[0].thermal_fuel_radius, sig_dig_fuel)

        beta = random.uniform(0, 1.0 + beta_delta)
        child_obj.cad_u_ratio         =  round(beta * individual_1[0].cad_u_ratio + (1 - beta) * individual_2[0].cad_u_ratio, sig_dig_plate)

        beta = random.uniform(0, 1.0 + beta_delta)
        child_obj.void_value          =  round(beta * individual_1[0].void_value + (1 - beta) * individual_2[0].void_value, sig_dig_plate)


        while child_obj.fast_fuel_radius > (maximum_fuel_radius):
            beta = random.uniform(0, 1.0 + beta_delta)
            child_obj.fast_fuel_radius = round(beta * individual_1[0].fast_fuel_radius + (1 - beta) * individual_2[0].fast_fuel_radius, sig_dig_fuel)
        while child_obj.thermal_fuel_radius > (maximum_fuel_radius):
            beta = random.uniform(0, 1.0 + beta_delta)
            child_obj.thermal_fuel_radius = round(beta * individual_1[0].thermal_fuel_radius + (1 - beta) * individual_2[0].thermal_fuel_radius, sig_dig_fuel)
        while child_obj.cad_u_ratio > 0.999 or child_obj.cad_u_ratio < 0.001:
            beta = random.uniform(0, 1.0 + beta_delta)
            child_obj.cad_u_ratio = round(beta * individual_1[0].cad_u_ratio + (1 - beta) * individual_2[0].cad_u_ratio, sig_dig_plate)
        while child_obj.void_value > 0.999 or child_obj.void_value < 0.001:
            beta = random.uniform(0, 1.0 + beta_delta)
            child_obj.void_value = round(beta * individual_1[0].void_value + (1 - beta) * individual_2[0].void_value, sig_dig_plate)

        return child_obj, children_list


    if crossover_type == 'averaging_crossover':

        print("Trying to crossover")
        zero_avg  = (individual_1[0].fast_fuel_radius + individual_2[0].fast_fuel_radius) / 2
        one_avg   = (individual_1[0].thermal_fuel_radius + individual_2[0].thermal_fuel_radius) / 2
        two_avg   = (individual_1[0].cad_u_ratio + individual_2[0].cad_u_ratio) / 2
        three_avg = (individual_1[0].void_value + individual_2[0].void_value) / 2





        if use_three_zones == True:
            four_avg = (individual_1[0].fast_fuel_2_radius + individual_2[0].fast_fuel_2_radius) / 2

        child_obj.fast_fuel_radius = round(zero_avg * random.uniform((1 - random_window), (1 + random_window)), sig_dig_fuel)
        child_obj.thermal_fuel_radius = round(one_avg * random.uniform((1 - random_window), (1 + random_window)), sig_dig_fuel)
        child_obj.cad_u_ratio = round(two_avg * random.uniform((1 - random_window), (1 + random_window)), sig_dig_plate)
        child_obj.void_value = round(three_avg * random.uniform((1 - random_window), (1 + random_window)), sig_dig_plate)
        if use_three_zones == True:
            child_obj.fast_fuel_2_radius = round(four_avg * random.uniform((1 - random_window), (1 + random_window)), sig_dig_fuel)
            while child_obj.fast_fuel_2_radius > (maximum_fuel_radius):
                child_obj.fast_fuel_2_radius = round(
                    four_avg * random.uniform((1 - random_window), (1 + random_window)),
                    4)

        while child_obj.fast_fuel_radius > (maximum_fuel_radius):
            child_obj.fast_fuel_radius = round(zero_avg * random.uniform((1 - random_window), (1 + random_window)), sig_dig_fuel)
        while child_obj.thermal_fuel_radius > (maximum_fuel_radius):
            child_obj.thermal_fuel_radius = round(one_avg * random.uniform((1 - random_window), (1 + random_window)), sig_dig_fuel)
        while child_obj.cad_u_ratio > 0.999:
            child_obj.cad_u_ratio = 0.999
        while child_obj.cad_u_ratio < 0.001:
            child_obj.cad_u_ratio = 0.001
        while child_obj.void_value > 0.999:
            child_obj.void_value = 0.999
        while child_obj.void_value < 0.001:
            child_obj.void_value = 0.001

        print("Parent 1 val: %s Parent 2 val: %s Child Val: %s" % (
            individual_1[0].fast_fuel_radius, individual_2[0].fast_fuel_radius, child_obj.fast_fuel_radius))

        # Job with just 1 fuel radius
        if job_type == 1:
            child_obj.thermal_fuel_radius = 0
            child_obj.cad_u_ratio = 0
            child_obj.void_value = 0
            child_obj.fast_fuel_2_radius = 0
        # Job with fast and thermal fuel radius
        if job_type == 2:
            child_obj.cad_u_ratio = 0
            child_obj.void_value = 0
            child_obj.fast_fuel_2_radius = 0
        # Job with fast and thermal fuel radius and plates
        if job_type == 3:
            child_obj.fast_fuel_2_radius = 0

        return child_obj, children_list

    if crossover_type == 'classic_crossover':
        # In this crossover method, the two parents will produce a child with a mix of the four variables of the parent
        # for example, if a random number of 75 is rolled, the resulting child will have 3 parts from parent 1,
        #  and 1 from parent 2
        print("Applying classic crossover technique")
        print("Parent 1 values: ", individual_1[0].fast_fuel_radius, individual_1[0].thermal_fuel_radius,
              individual_1[0].cad_u_ratio, individual_1[0].void_value)
        print("Parent 2 values: ", individual_2[0].fast_fuel_radius, individual_2[0].thermal_fuel_radius,
              individual_2[0].cad_u_ratio, individual_2[0].void_value)
        split_value = random.randint(1, 100)
        if split_value <= 25:
            child_obj.fast_fuel_radius = individual_1[0].fast_fuel_radius
            child_obj.thermal_fuel_radius = individual_2[0].thermal_fuel_radius
            child_obj.cad_u_ratio = individual_2[0].cad_u_ratio
            child_obj.void_value = individual_2[0].void_value

        if split_value <= 50:
            child_obj.fast_fuel_radius = individual_1[0].fast_fuel_radius
            child_obj.thermal_fuel_radius = individual_1[0].thermal_fuel_radius
            child_obj.cad_u_ratio = individual_2[0].cad_u_ratio
            child_obj.void_value = individual_2[0].void_value

        if split_value <= 75:
            child_obj.fast_fuel_radius = individual_1[0].fast_fuel_radius
            child_obj.thermal_fuel_radius = individual_1[0].thermal_fuel_radius
            child_obj.cad_u_ratio = individual_1[0].cad_u_ratio
            child_obj.void_value = individual_2[0].void_value

        if split_value <= 100:
            child_obj.fast_fuel_radius = individual_1[0].fast_fuel_radius
            child_obj.thermal_fuel_radius = individual_1[0].thermal_fuel_radius
            child_obj.cad_u_ratio = individual_1[0].cad_u_ratio
            child_obj.void_value = individual_2[0].void_value
        child_params_string = str(child_obj.fast_fuel_radius) + \
                              "," + str(child_obj.thermal_fuel_radius) + \
                              "," + str(child_obj.cad_u_ratio) + \
                              "," + str(child_obj.void_value)


    # Job with just 1 fuel radius
    if job_type == 1:
        child_obj.thermal_fuel_radius = 0
        child_obj.cad_u_ratio = 0.5
        child_obj.void_value = 1.0
        child_obj.fast_fuel_2_radius = 0
    # Job with fast and thermal fuel radius
    if job_type == 2:
        child_obj.cad_u_ratio = 0.5
        child_obj.void_value = 1.0
        child_obj.fast_fuel_2_radius = 0
    # Job with fast and thermal fuel radius and plates
    if job_type == 3:
        child_obj.fast_fuel_2_radius = 0
    print(child_params_string)
    for item in children_list:
        if item == child_params_string:
            print("FOUND A COPY!")
        while item == child_params_string:
            print("This child has already been created. Mutating it. ")
            split_value_2 = random.randint(1, 100)
            if split_value_2 <= 25:
                # mutating the value based on random_window specification
                child_obj.fast_fuel_radius = round(child_obj.fast_fuel_radius *
                                                   random.uniform((1 - random_window),
                                                                  (1 + random_window)), sig_dig_fuel)
                # testing to make sure the value is within the specifications of the model
                while child_obj.fast_fuel_radius > (maximum_fuel_radius):
                    child_obj.fast_fuel_radius = round(
                        child_obj.fast_fuel_radius * random.uniform((1 - random_window), (1 + random_window)), sig_dig_fuel)
            if split_value_2 > 25 and split_value_2 <= 50:
                child_obj.thermal_fuel_radius = round(child_obj.thermal_fuel_radius *
                                                      random.uniform((1 - random_window),
                                                                     (1 + random_window)), sig_dig_fuel)
                while child_obj.thermal_fuel_radius > (maximum_fuel_radius):
                    child_obj.thermal_fuel_radius = round(
                        child_obj.thermal_fuel_radius * random.uniform((1 - random_window), (1 + random_window)), sig_dig_fuel)
            if split_value_2 > 50 and split_value_2 <= 75:
                child_obj.cad_u_ratio = round(child_obj.cad_u_ratio * random.uniform((1 - random_window),
                                                                                     (1 + random_window)), sig_dig_plate)
                while child_obj.cad_u_ratio > 0.999:
                    child_obj.cad_u_ratio = 0.999
                while child_obj.cad_u_ratio < 0.001:
                    child_obj.cad_u_ratio = 0.001
            if split_value_2 > 75 and split_value_2 <= 100:
                child_obj.void_value = round(child_obj.void_value * random.uniform((1 - random_window),
                                                                                   (1 + random_window)), sig_dig_plate)
                while child_obj.void_value > 0.999:
                    child_obj.void_value = 0.999
                while child_obj.void_value < 0.001:
                    child_obj.void_value = 0.001
            child_params_string = str(child_obj.fast_fuel_radius) + \
                                  "," + str(child_obj.thermal_fuel_radius) + \
                                  "," + str(child_obj.cad_u_ratio) + \
                                  "," + str(child_obj.void_value)
            print("new: ", child_params_string)
    if child_params_string not in children_list:
        children_list.append(child_params_string)
    print("Split value: ", split_value)
    print("Child values: ", child_obj.fast_fuel_radius, child_obj.thermal_fuel_radius,
          child_obj.cad_u_ratio, child_obj.void_value)

    return child_obj, children_list


toolbox.register("crossover", crossover)


def write_from_restart_file(restart_file_string, pop):
    print("Overwriting initial jobs with info from restart file")
    restart_count_int = 1
    line_list = []
    with open(restart_file_string, 'r') as restart_file:
        for line in restart_file:
            if line.startswith('c'): continue
            line_list.append(line.strip())
    for restart in line_list:
        restart_split = restart.split(',')
        print(restart_split)
        pop[int(restart_split[0])][0].fast_fuel_radius = float(restart_split[1])
        pop[int(restart_split[0])][0].thermal_fuel_radius = float(restart_split[2])
        pop[int(restart_split[0])][0].cad_u_ratio = float(restart_split[3])
        pop[int(restart_split[0])][0].void_value = float(restart_split[4])
    return pop


def waiting_on_cluster_jobs(list_of_jobs_waiting_on, job_type, completed_run_tag):
    if list_of_jobs_waiting_on == []:
        print("No jobs to wait on, continuing.")
        return
    continue_flag = True
    gen_str = "gen_" + str(g) + "_"
    job_count = 0

    for job in list_of_jobs_waiting_on:
        job_count += 1
    print("Waiting on ", job_count, " jobs")
    while continue_flag == True:
        finished_job = 0
        list_of_files = os.listdir()
        for file in list_of_files:
            if gen_str not in file:
                continue
            if job_type not in file:
                continue
            if completed_run_tag not in file:
                continue
            finished_job += 1
            if finished_job == job_count:
                continue_flag = False
        if continue_flag == True:
            print("Waiting 15 seconds", "jobs finished:", finished_job, "of", job_count)
            time.sleep(15)

def bin_individuals_by_keff(pop):
    print("Collecting jobs for source run. A minimum of ", str(number_of_source_jobs_to_run), "will be run.")
    list_of_jobs = []
    keff_max_local = keff_max
    keff_min_local = keff_min
    number_of_source_jobs_found = 0
    continue_to_gather_source_jobs = True
    gen_str = "_gen_" + str(generation)
    #Gathering keffs for each job
    for ind in pop:
        input_filename = ind[0].input_name
        if run_on_necluster == True:
            keff = get_keff(input_filename + "_keff.inpo")
            ind[0].keff = keff
        if run_on_necluster == False:
            keff = get_keff(input_filename + "_keff.inp.out")
            ind[0].keff = keff

    while continue_to_gather_source_jobs:
        for ind in pop:
            keff = ind[0].keff
            if float(keff) <= keff_max_local:
                if float(keff) > keff_min_local:
                    if ind[0].ran_true_source == True:
                        continue
                    if gen_str not in ind[0].input_name:
                        continue
                    list_of_jobs.append(ind[0].input_name)
                    ind[0].ran_true_source = True
                    number_of_source_jobs_found += 1
        if number_of_source_jobs_found > number_of_source_jobs_to_run:
            continue_to_gather_source_jobs = False

        if continue_to_gather_source_jobs == True:
            keff_max_local = keff_max_local + 0.0005
            keff_min_local = keff_min_local - 0.0005

            if keff_max_local >= maximum_keff_of_source_job:
                print("Maximum range keff of",str(maximum_keff_of_source_job),"has been reached. Continuing with",str(number_of_source_jobs_found),"jobs")
                return pop, list_of_jobs

            print("Found",str(number_of_source_jobs_found), "jobs. Expanding range to", str(keff_max_local), str(keff_min_local))
            continue
        print("Found ", str(number_of_source_jobs_found), "jobs. In the bin:",str(keff_max_local), str(keff_min_local))
    return pop, list_of_jobs


g = 0
# fits_keff = [ind[0].fitness.values[0] for ind in pop]
# fits_mass = [ind[0].fitness.values[1] for ind in pop]

print("Creating output file")
output_file = open(output_file_str, 'w')
output_file.write(
    "ind.input_filename,ind.keff,ind.mass,ind.fast_fuel_radius,ind.fast_fuel_2_radius,ind.thermal_fuel_radius,ind.cad_u_ratio,ind.void_value,fast_flux,thermal_flux,ind.fitness.values[0],ind.fitness.values[1],ind.fitness.values[2],ind.fitness.values[3]\n")
output_file.close()

print("Building " + str(population_int) + " initial models.")
pop = toolbox.population(population_int)

if restart_from_file == True:
    pop = write_from_restart_file(restart_file_string, pop)

# Building initial mcnp models
list(map(toolbox.build_mcnp_model, pop))
# Running keff calc on initial models
print("Running the initial keff jobs.")

if run_on_necluster == False:
    toolbox.run_Mcnp("keff", g, "keff")
if run_on_necluster == True:
    run_file = run_on_cluster("keff", g, "keff")
    print(run_file)
    waiting_on_cluster_jobs(run_file, "keff", "_done.dat")

print("Running true source calc if %s < keff < %s" % (str(keff_min), str(keff_max)))
true_source_run_list = []

pop, true_source_run_list = bin_individuals_by_keff(pop)

#for ind in pop:
#    input_filename = ind[0].input_name
#    if run_on_necluster == True:
#        keff = get_keff(input_filename + "_keff.inpo")
#    if run_on_necluster == False:
#        keff = get_keff(input_filename + "_keff.inp.out")
#    ind[0].keff = keff
#    if float(keff) <= keff_max:
#        if float(keff) > keff_min:
#            print("Job %s has a keff of %s, so running true source calc" % (str(ind[0].input_name), str(ind[0].keff)))
#            # Setting "Ran true source" to true for this individual
#            ind[0].ran_true_source = True
#            true_source_run_list.append(ind[0].input_name)
#
#            # running the true source calcs


for job in true_source_run_list:
    if run_on_necluster == True:
        run_on_cluster("source", g, job)
    if run_on_necluster == False:
        toolbox.run_Mcnp("source", g, job)
waiting_on_cluster_jobs(true_source_run_list, "source", "_done.dat")

while g < number_of_generations:
    output_file = open(output_file_str, 'a')
    # Updating local generation variable
    g = g + 1
    # Updating global generation variable
    generation = generation + 1
    print("-- Generation %i --" % g)

    print("Evaluating Last Gen Jobs")
    # evaluating initial population, keff, mass, and flux
    fitnesses = list(map(toolbox.evaluate, pop))
    for ind, fit in zip(pop, fitnesses):
        print("ind, fit:", ind[1], fit)
        ind.fitness.values = fit

    print("Best Guess:", pop[0][0].fast_fuel_radius, pop[0][0].thermal_fuel_radius, pop[0][0].cad_u_ratio,
          pop[0].fitness)

    pop = sorted(pop, key=lambda pop: pop.fitness, reverse=True)

    # Writing out job and score for previous generation to file
    total_keff = 0.0
    total_mass = 0.0
    for ind in pop:
        output_file.write(str(ind[0].input_name) + "," + str(ind[0].keff) + "," + str(ind[0].mass) + "," + str(
            ind[0].fast_fuel_radius) + "," + str(ind[0].fast_fuel_2_radius) + "," + str(
            ind[0].thermal_fuel_radius) + "," + str(
            ind[0].cad_u_ratio) + "," + str(ind[0].void_value) + "," + str(ind[0].fast_flux) + "," + str(
            ind[0].thermal_flux) + "," +
                          str(ind.fitness.values[0]) + "," + str(ind.fitness.values[1]) + "," + str(
            ind.fitness.values[2]) + ","
                          + str(ind.fitness.values[3]) + "\n")
        total_keff = total_keff + float(ind[0].keff)
        total_mass = total_mass + float(ind[0].mass)
    output_file.close()
    average_keff = total_keff / float(population_int)
    average_mass = total_mass / float(population_int)

    # Selecting the next generation individuals
    print("Selecting the parents of next generation")
    parents = toolbox.select(pop, keepers_int)

    # adding the selected files to the pass list, which won't be deleted in the cleanup function
    pass_list = []
    for ind in parents:
        pass_list.append(ind[1])

    # Clone the selected individuals
    print("Clone the selected individuals???")
    # parents = list(map(toolbox.clone, parents))
    # print(parents)

    # Apply crossover and mutation on the offspring
    print("Apply crossover and mutation on the offspring")
    count = 0
    parent_int_max = keepers_int - 1
    print("The following are the original fast and thermal diams")
    # for ind in pop:
    # print(ind[0].fast_fuel_radius, ind[0].thermal_fuel_radius)
    # Applying crossover

    for individual in pop:
        # getting int for parent
        parent_1 = random.randint(0, parent_int_max)
        parent_2 = random.randint(0, parent_int_max)

        # verifying 2 different parents used
        while parent_1 == parent_2:
            parent_2 = random.randint(0, parent_int_max)

        # keeping the top # from previous run
        if count <= parent_int_max:
            pop[count] = parents[count]
            pop[count].is_parent = True
            count += 1
            continue

        # combining parents into new ind
        print("Applying crossover to %s" % str(count))
        pop[count][0], children_list = crossover(parents[parent_1], parents[parent_2], children_list)
        pop[count][1] = job_name()
        del individual.fitness.values
        count += 1
    print(children_list)
    # for ind in pop:
    #   print(ind[0].fast_fuel_radius, ind[0].thermal_fuel_radius)

    # Chance to mutate the offspring
    print("Chance to mutate the offspring")
    count = 0
    for mutant in pop:
        # Not going to mutate the parents
        if count < keepers_int:
            count += 1
            continue
        if random.uniform(0.0, 1.0) < MUTPB:
            print("Mutation in ", mutant[1])
            toolbox.mutate(mutant)
            del mutant.fitness.values
        count += 1
    # checking for duplicate children
    pop = check_duplicate(pop)
    # Building new mcnp models
    print("Building mcnp models for this gen")
    list(map(toolbox.build_mcnp_model, pop))

    # Running mcnp on jobs
    print("Running mcnp keff jobs")
    if run_on_necluster == True:
        run_file = run_on_cluster("keff", g, "keff")
        print(run_file)
        waiting_on_cluster_jobs(run_file, "keff", "_done.dat")

    if run_on_necluster == False:
        toolbox.run_Mcnp("keff", g, "keff")
    print("Running true source calc if %s < keff < %s" % (str(keff_min), str(keff_max)))
    true_source_run_list = []

    #This function bins the jobs by keff
    pop, true_source_run_list = bin_individuals_by_keff(pop)


    # running the true source calcs
    for job in true_source_run_list:
        if run_on_necluster == True:
            run_on_cluster("source", g, job)

        if run_on_necluster == False:
            toolbox.run_Mcnp("source", g, job)
    waiting_on_cluster_jobs(true_source_run_list, "source", "_done.dat")