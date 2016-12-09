\* Cutting stock Master *\
Minimize
OBJ: 50 occurrence_occ_0 + 50 occurrence_occ_1 + 50 occurrence_occ_2
 + 100 occurrence_occ_3 + 100 occurrence_occ_4 + 100 occurrence_occ_5
Subject To
demand_for_order_0_should_be_met: occurrence_occ_0 + 2 occurrence_occ_3
 + 2 occurrence_occ_4 + 2 occurrence_occ_5 = 4
demand_for_order_1_should_be_met: occurrence_occ_1 = 2
demand_for_order_2_should_be_met: occurrence_occ_2 = 0
demand_for_order_3_should_be_met: 5 occurrence_occ_3 + occurrence_occ_5 = 2
demand_for_order_4_should_be_met: 5 occurrence_occ_4 = 3
demand_for_order_5_should_be_met: 2 occurrence_occ_5 = 1
End
