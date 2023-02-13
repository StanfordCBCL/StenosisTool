from svinterface.core.zerod.lpn import OriginalLPN, LPN

if __name__ == '__main__':
    lpn = OriginalLPN.from_file("test.in")

    lpn = OriginalLPN.from_dict(lpn.lpn_data)
    test1 = lpn.get_vessel("branch0_seg0")
    lpn.change_vessel("branch0_seg0", R = 10, C = 11, L = 12, S = 13)
    test3 = lpn.get_vessel(0)
    
    #print(lpn.lpn_data)
    
    
    test1 = lpn.get_junction(0)
    test2 = lpn.get_junction("J0")
    
    #print(test1, test2)
    
    fast = lpn.get_fast_lpn()
    
    #print(fast.data)
    
    slow = lpn.get_full_lpn()
    
    print(slow.lpn_data.keys())
    
    print(slow.bc_data.keys())

    tmp = """LPA
LPA_1
LPA_10
LPA_11
LPA_12
LPA_13
LPA_14
LPA_15
LPA_16
LPA_17
LPA_18
LPA_19
LPA_2
LPA_20
LPA_21
LPA_22
LPA_23
LPA_24
LPA_25
LPA_26
LPA_27
LPA_28
LPA_29
LPA_3
LPA_30
LPA_31
LPA_32
LPA_33
LPA_34
LPA_35
LPA_36
LPA_37
LPA_38
LPA_39
LPA_4
LPA_40
LPA_41
LPA_42
LPA_43
LPA_44
LPA_45
LPA_46
LPA_5
LPA_6
LPA_7
LPA_8
LPA_9
RPA
RPA_1
RPA_10
RPA_11
RPA_12
RPA_13
RPA_14
RPA_15
RPA_16
RPA_17
RPA_18
RPA_19
RPA_2
RPA_20
RPA_21
RPA_22
RPA_23
RPA_24
RPA_25
RPA_26
RPA_27
RPA_28
RPA_29
RPA_3
RPA_30
RPA_31
RPA_32
RPA_33
RPA_34
RPA_35
RPA_36
RPA_37
RPA_38
RPA_39
RPA_4
RPA_40
RPA_41
RPA_42
RPA_43
RPA_44
RPA_45
RPA_46
RPA_5
RPA_6
RPA_7
RPA_8
RPA_9
""".split("\n")

    slow.add_rcrt_map(tmp)
    
    print(slow.get_outlet_bc("RCR_0"))
    print(slow.flags)
    
    slow.write_lpn_file("test2.out")

    new_slow = LPN.from_file("test2.out")
    print(new_slow.get_outlet_bc("RCR_0"))
    print(new_slow.get_vessel(0))
    mpa = slow.get_tree()
    
    counter = 0
    for node in slow.tree_bfs_iterator(mpa):
        print(node)
        counter+=1
    print("Total:", counter)
    
    counter = 0
    for node in slow.tree_bfs_iterator(mpa, "branch"):
        print(node)
        counter+=1
    print("Total:", counter)
    
    counter = 0
    for node in slow.tree_bfs_iterator(mpa, "junction"):
        print(node)
        counter+=1
    print("Total:", counter)
        
        
    slow.det_lpa_rpa(mpa)
    slow.write_lpn_file("test3.out")
    for node in slow.tree_bfs_iterator(mpa):
        for vess in node.vessel_info:
            print(vess['side'])
    

    new_slow = LPN.from_file("test3.out")
    new_mpa = new_slow.get_tree()
    print(new_mpa.vessel_info[0]['side'])
    
    
