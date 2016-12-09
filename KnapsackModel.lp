\* Knapsack problem *\
Maximize
OBJ: 50 assign_a_0 + 50 assign_a_1
Subject To
do_not_fill_more_than_capacity_for_geo_0: 25 assign_a_0 + 25 assign_a_1
 + 25 assign_a_2 <= 5
do_not_fill_more_than_capacity_for_geo_1: 10 assign_a_3 + 10 assign_a_4
 + 20 assign_a_5 <= 5
Bounds
0 <= assign_a_0
0 <= assign_a_1
0 <= assign_a_2
0 <= assign_a_3
0 <= assign_a_4
0 <= assign_a_5
Generals
assign_a_0
assign_a_1
assign_a_2
assign_a_3
assign_a_4
assign_a_5
End
