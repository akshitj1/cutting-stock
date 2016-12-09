import copy, json, time
from ColGen import col_gen
import numpy as np

import sys
sys.path.append("/home/kiran/Adplanner/ad_planner/schedule_planner_python")
from generate_day_part_plan import generate_dat_part_plan


def initial_soln(ad_planner_obj, input_data):

    output = ad_planner_obj.max_margin_branch_across_processes[0]
    geo_list = input_data.keys()
    geo_list.sort()
    if "Master" in geo_list:
        geo_list.remove("Master")
        geo_list.insert(0, "Master")

    duration, min_rotates, mcc, captions = [], [], [], []
    geo_idx = ()
    st = 0
    for geo in geo_list:
        duration += input_data[geo]['duration']
        min_rotates += input_data[geo]['min_rotates']
        mcc += input_data[geo]['multiple_caption_combination']
        captions += input_data[geo]['captions']
        en = st + len(input_data[geo]['captions']) - 1
        geo_idx += ((st, en),)
        st = en + 1

    B = np.identity(len(duration))
    cost = np.zeros(len(duration))

    col_idx = 0
    for branch in output:
        pattern = []
        for geo_pattern in branch:
            pattern += list(geo_pattern)

        B[:, col_idx] = np.array(pattern)
        cost[col_idx] = ad_planner_obj.best_spot_length[col_idx]
        col_idx += 1

    cost[col_idx:] = duration[col_idx:]
    cost = np.matrix(cost)

    return duration, min_rotates, mcc, captions, geo_idx, B, cost


def ensure_same_geos_across_dayparts(all_day_part_data):
    geos = set([])
    example_daypart_object = {}
    # To add all geos in the input to the set geos.
    for day_part, day_part_data in all_day_part_data.iteritems():
        geos.update(day_part_data['main_data'].keys())
        # Getting an example/sample object to create an empty object.
        if not example_daypart_object:
            if day_part_data['main_data']:
                example_daypart_object = copy.deepcopy(day_part_data['main_data'][day_part_data['main_data'].keys()[0]])

    if example_daypart_object:
        for key in example_daypart_object:
            example_daypart_object[key] = []

        for day_part in all_day_part_data:
            for geo in geos - set(all_day_part_data[day_part]['main_data'].keys()):
                all_day_part_data[day_part]['main_data'][geo] = example_daypart_object

    return all_day_part_data


def getInputFromPlan(filename="plan-1.json", selectedDayparts = []):

    text = open(filename).read()
    allData = json.loads(text)
    if 'algorithmInput' in allData.keys():
        obj = allData['algorithmInput']
    elif 'generatedInput' in allData.keys():
        obj = allData['generatedInput']
    else:
        obj = allData
    obj['job_id'] = str(int(time.time()))+"-DEVELOPPER"
    text = json.dumps(obj)
    return text


def get_input_args(input_data):

    all_geos = input_data.keys()[0:4]
    duration, min_rotates, mcc, captions = [], [], [], []
    geo_idx = ()
    st = 0
    for geo in all_geos:
        duration += input_data[geo]['duration']
        min_rotates += input_data[geo]['min_rotates']
        mcc += input_data[geo]['multiple_caption_combination']
        captions += input_data[geo]['captions']
        en = st + len(input_data[geo]['captions']) - 1
        geo_idx += ((st, en),)
        st = en + 1

    return duration, min_rotates, mcc, captions, geo_idx

inp = getInputFromPlan("/home/kiran/Scenario Builder/Input/Algo-2-Suboptimal/ZT_2016_05_01.json", "")

all_data = json.loads(inp)
all_day_part_data = all_data['all_day_part_data']
input_data = all_day_part_data["Morning"]

input_data['path'] = '/home/kiran/Adplanner/ad_planner'


start_time = time.time()
ad_planner_obj = generate_dat_part_plan(input_data)

input_data = input_data['main_data']

duration, min_rotates, mcc, captions, geo_idx, B, cost = initial_soln(ad_planner_obj, input_data)

B, C, Xb = col_gen(duration, min_rotates, geo_idx, B, cost)
print np.round(Xb, 4)