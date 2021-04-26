function [names,T_c_arr] = ising_names(n)
    switch n
        case 1
            
            names = {{
                'Ising2D_g=0.0000e+00_chi=8_02_April_2021_11:46';
                'Ising2D_g=0.0000e+00_chi=11_02_April_2021_11:46';
                'Ising2D_g=0.0000e+00_chi=16_02_April_2021_11:46';
                'Ising2D_g=0.0000e+00_chi=23_02_April_2021_11:46';
                'Ising2D_g=0.0000e+00_chi=32_02_April_2021_11:46';
                'Ising2D_g=0.0000e+00_chi=45_02_April_2021_11:46';
                },{
                'Ising2D_g=1.5000e+00_chi=8_02_April_2021_11:17';
                'Ising2D_g=1.5000e+00_chi=11_02_April_2021_11:17';
                'Ising2D_g=1.5000e+00_chi=16_02_April_2021_11:17';
                'Ising2D_g=1.5000e+00_chi=23_02_April_2021_11:17';
                'Ising2D_g=1.5000e+00_chi=32_02_April_2021_11:17';
                'Ising2D_g=1.5000e+00_chi=45_02_April_2021_11:17';
                'Ising2D_g=1.5000e+00_chi=64_02_April_2021_11:17';
                },{
                'Ising2D_g=2.5000e+00_chi=8_02_April_2021_09:38';
                'Ising2D_g=2.5000e+00_chi=11_02_April_2021_09:38';
                'Ising2D_g=2.5000e+00_chi=16_02_April_2021_09:38';
                'Ising2D_g=2.5000e+00_chi=23_02_April_2021_09:38';
                'Ising2D_g=2.5000e+00_chi=32_02_April_2021_09:38';
                'Ising2D_g=2.5000e+00_chi=45_02_April_2021_09:38';
                },{
                'Ising2D_g=2.5000e+00_chi=20_22_March_2021_10:11_a';
                'Ising2D_g=2.5000e+00_chi=25_22_March_2021_10:11_a';
                'Ising2D_g=2.5000e+00_chi=30_22_March_2021_10:11_a';
                'Ising2D_g=2.5000e+00_chi=40_22_March_2021_10:46_a';
                'Ising2D_g=2.5000e+00_chi=60_22_March_2021_11:49_a';
                'Ising2D_g=2.5000e+00_chi=65_22_March_2021_13:52_a';
                'Ising2D_g=2.5000e+00_chi=70_22_March_2021_13:51_a';
                }};
            
            T_c_arr= [2.2691859, 1.97951, 1.27372, 0.9];

        case 2

            names = {{
                'Ising2D_g=0.0000e+00_chi=8_08_April_2021_09:52';
                'Ising2D_g=0.0000e+00_chi=11_08_April_2021_10:11';
                'Ising2D_g=0.0000e+00_chi=16_08_April_2021_10:43';
                'Ising2D_g=0.0000e+00_chi=23_08_April_2021_11:50';
                'Ising2D_g=0.0000e+00_chi=32_08_April_2021_14:41';
                'Ising2D_g=0.0000e+00_chi=45_08_April_2021_18:59';
                'Ising2D_g=0.0000e+00_chi=64_09_April_2021_05:59';
                'Ising2D_g=0.0000e+00_chi=91_11_April_2021_18:25';
                'Ising2D_g=0.0000e+00_chi=128_16_April_2021_16:04';
                }, {
                'Ising2D_g=1.5000e+00_chi=8_08_April_2021_10:15';
                'Ising2D_g=1.5000e+00_chi=11_08_April_2021_10:28';
                'Ising2D_g=1.5000e+00_chi=16_08_April_2021_10:48';
                'Ising2D_g=1.5000e+00_chi=23_08_April_2021_11:29';
                'Ising2D_g=1.5000e+00_chi=32_08_April_2021_12:19';
                'Ising2D_g=1.5000e+00_chi=45_08_April_2021_16:50';
                'Ising2D_g=1.5000e+00_chi=64_09_April_2021_04:35';
                'Ising2D_g=1.5000e+00_chi=91_10_April_2021_20:23';
                'Ising2D_g=1.5000e+00_chi=128_13_April_2021_18:47';
                }, {
                'Ising2D_g=2.5000e+00_chi=8_08_April_2021_09:52';
                'Ising2D_g=2.5000e+00_chi=11_08_April_2021_10:02';
                'Ising2D_g=2.5000e+00_chi=16_08_April_2021_10:17';
                'Ising2D_g=2.5000e+00_chi=23_08_April_2021_10:52';
                'Ising2D_g=2.5000e+00_chi=32_08_April_2021_12:39';
                'Ising2D_g=2.5000e+00_chi=45_08_April_2021_17:46';
                'Ising2D_g=2.5000e+00_chi=64_09_April_2021_05:52';
                'Ising2D_g=2.5000e+00_chi=91_10_April_2021_23:23';
                'Ising2D_g=2.5000e+00_chi=128_14_April_2021_21:21';
                }, {
                'Ising2D_g=2.9000e+00_chi=8_08_April_2021_09:52';
                'Ising2D_g=2.9000e+00_chi=11_08_April_2021_10:03';
                'Ising2D_g=2.9000e+00_chi=16_08_April_2021_10:06';
            }};
         T_c_arr = [2.2691859, 1.97951, 1.27372, 0.9];
    end

end