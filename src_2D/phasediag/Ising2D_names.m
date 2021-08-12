function [names, T_c_arr] = Ising2D_names(n)
    switch n
        case 1

            names = {{
                'Ising2D_g=0.0000e+00_chi=8_02_April_2021_11:46';
                'Ising2D_g=0.0000e+00_chi=11_02_April_2021_11:46';
                'Ising2D_g=0.0000e+00_chi=16_02_April_2021_11:46';
                'Ising2D_g=0.0000e+00_chi=23_02_April_2021_11:46';
                'Ising2D_g=0.0000e+00_chi=32_02_April_2021_11:46';
                'Ising2D_g=0.0000e+00_chi=45_02_April_2021_11:46';
                }, {
                'Ising2D_g=1.5000e+00_chi=8_02_April_2021_11:17';
                'Ising2D_g=1.5000e+00_chi=11_02_April_2021_11:17';
                'Ising2D_g=1.5000e+00_chi=16_02_April_2021_11:17';
                'Ising2D_g=1.5000e+00_chi=23_02_April_2021_11:17';
                'Ising2D_g=1.5000e+00_chi=32_02_April_2021_11:17';
                'Ising2D_g=1.5000e+00_chi=45_02_April_2021_11:17';
                'Ising2D_g=1.5000e+00_chi=64_02_April_2021_11:17';
                }, {
                'Ising2D_g=2.5000e+00_chi=8_02_April_2021_09:38';
                'Ising2D_g=2.5000e+00_chi=11_02_April_2021_09:38';
                'Ising2D_g=2.5000e+00_chi=16_02_April_2021_09:38';
                'Ising2D_g=2.5000e+00_chi=23_02_April_2021_09:38';
                'Ising2D_g=2.5000e+00_chi=32_02_April_2021_09:38';
                'Ising2D_g=2.5000e+00_chi=45_02_April_2021_09:38';
                }, {
                'Ising2D_g=2.5000e+00_chi=20_22_March_2021_10:11_a';
                'Ising2D_g=2.5000e+00_chi=25_22_March_2021_10:11_a';
                'Ising2D_g=2.5000e+00_chi=30_22_March_2021_10:11_a';
                'Ising2D_g=2.5000e+00_chi=40_22_March_2021_10:46_a';
                'Ising2D_g=2.5000e+00_chi=60_22_March_2021_11:49_a';
                'Ising2D_g=2.5000e+00_chi=65_22_March_2021_13:52_a';
                'Ising2D_g=2.5000e+00_chi=70_22_March_2021_13:51_a';
                }};

            T_c_arr = [2.2691859, 1.97951, 1.27372, 0.9];

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
            %'TIM_g=2.5_order_7_chi=8_trunc_20_sym=1_05_June_2021_09:24';
                }};
            T_c_arr = [2.2691859, 1.97951, 1.27372];

        case 3

            names = {{
                'TIM_T=0.7_order_5_chi=8_trunc_20_sym=1_21_May_2021_15:53';
                'TIM_T=0.7_order_5_chi=11_trunc_20_sym=1_21_May_2021_15:57';
                'TIM_T=0.7_order_5_chi=16_trunc_20_sym=1_21_May_2021_16:26';
                'TIM_T=0.7_order_5_chi=23_trunc_20_sym=1_21_May_2021_16:28';
            %'TIM_T=0.7_order_6_chi=8_trunc_20_sym=1_22_May_2021_10:09'; %L=6 bd=20
            %'TIM_T=0.7_order_5_chi=8_trunc_20_sym=1_25_May_2021_17:35';
            %% complex version
                }};
            T_c_arr = [2.817];

        case 4
            names = {{
                'TIM_T=0.7_order_6_chi=8_trunc_20_sym=1_22_May_2021_10:09'; %L=6 bd=20
                'TIM_T=0.7_order_6_chi=11_trunc_20_sym=1_22_May_2021_14:45';
                'TIM_T=0.7_order_6_chi=16_trunc_20_sym=1_22_May_2021_14:45';
                'TIM_T=0.7_order_6_chi=23_trunc_20_sym=1_22_May_2021_20:26';
                'TIM_T=0.7_order_6_chi=8_trunc_30_sym=1_03_June_2021_08:37';
                }};

            T_c_arr = [2.8544];

        case 5
            names = {{
            %'TIM_g=2.9_order_6_chi=6_trunc_30_sym=1_03_June_2021_08:33';
            %'TIM_g=2.9_order_6_chi=6_trunc_20_sym=1_03_June_2021_09:59'; %no loops
                'TIM_g=2.9_order_7_chi=6_trunc_64_sym=1_03_June_2021_09:41';
                'TIM_g=2.9_order_7_chi=11_trunc_64_sym=1_03_June_2021_11:54';
                'TIM_g=2.9_order_7_chi=16_trunc_64_sym=1_03_June_2021_11:56';
                'TIM_g=2.9_order_7_chi=23_trunc_64_sym=1_03_June_2021_20:14';
            %'TIM_g=2.9_order_7_chi=32_trunc_64_sym=1_04_June_2021_08:50';
            %'TIM_g=2.9_order_7_chi=8_trunc_30_sym=1_05_June_2021_10:03';
            %'TIM_g=2.9_order_7_chi=11_trunc_30_sym=1_05_June_2021_11:08';
                }};

            T_c_arr = [0.605];
        case 6
            names = {{
            %'TIM_g=2.5_order_5_chi=8_sym=1_30_July_2021_18:26';
                'TIM_g=2.5_order_5_chi=11_sym=1_31_July_2021_10:08';
                'TIM_g=2.5_order_5_chi=16_sym=1_31_July_2021_10:08';
                'TIM_g=2.5_order_5_chi=23_sym=1_31_July_2021_12:14';
                'TIM_g=2.5_order_5_chi=32_sym=1_31_July_2021_14:56';
                'TIM_g=2.5_order_5_chi=45_sym=1_01_August_2021_14:01';
                'TIM_g=2.5_order_5_chi=64_sym=1_01_August_2021_14:00';
                }};
            T_c_arr = [1.2736];
        case 7 %6 but reprocessed
            names = {{
                'TIM_g=2.5_order_5_chi=11_sym=1_31_July_2021_10:08_1x1';
                'TIM_g=2.5_order_5_chi=16_sym=1_31_July_2021_10:08_1x1';
                'TIM_g=2.5_order_5_chi=23_sym=1_31_July_2021_12:14_1x1';
                'TIM_g=2.5_order_5_chi=32_sym=1_31_July_2021_14:56_1x1';
                'TIM_g=2.5_order_5_chi=45_sym=1_01_August_2021_14:01_1x1';
                'TIM_g=2.5_order_5_chi=64_sym=1_01_August_2021_14:00_1x1';
                }};
            T_c_arr = [1.2736];
    end

end
