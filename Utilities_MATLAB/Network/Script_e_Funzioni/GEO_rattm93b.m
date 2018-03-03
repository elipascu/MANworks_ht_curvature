% 		
% nodes=[1	200	250	100	0	0.4056	;
% 2	170	215	100	2	0	;
% 3	115	190	100	2	0	;
% 4	60	165	100	2	0	;
% 5	95	153	140	2	0	;
% 6	125	142	140	2	0	;
% 7	100	65	100	2	0	;
% 8	85	0	100	0	0	;
% 9	125	55	100	2	0	;
% 10	160	80	100	2	0	;
% 11	205	70	100	2	0	;
% 12	160	130	100	2	0	;
% 13	370	40	100	0	0.4056	;
% 14	255	0	100	0	0	;
% 15	240	170	100	2	0	;
% 16	270	145	100	2	0	;
% 17	300	185	100	2	0	;
% 18	280	215	100	2	0	;
% 19	300	250	100	0	0.4056	;
% 20	265	70	100	2	0	;
% 21	315	80	100	2	0	;
% 22	330	125	100	2	0	;
% 23	370	135	100	0	0.4056	;
% 24	0	155	100	0	0	;
% 25	0	100	100	0	0	;
% 26	340	250	100	0	0.4056]	;

nodes=[1	4	5	2	0	0.4056	;
2	3.4	4.3	2	2	0	;
3	2.3	3.8	2	2	0	;
4	1.2	3.3	2	2	0	;
5	1.9	3.06	2.8	2	0	;
6	2.5	2.84	2.8	2	0	;
7	2	1.3	2	2	0	;
8	1.7	0	2	0	0	;
9	2.5	1.1	2	2	0	;
10	3.2	1.6	2	2	0	;
11	4.1	1.4	2	2	0	;
12	3.2	2.6	2	2	0	;
13	7.4	0.8	2	0	0.4056	;
14	5.1	0	2	0	0	;
15	4.8	3.4	2	2	0	;
16	5.4	2.9	2	2	0	;
17	6	3.7	2	2	0	;
18	5.6	4.3	2	2	0	;
19	6	5	2	0	0.4056	;
20	5.3	1.4	2	2	0	;
21	6.3	1.6	2	2	0	;
22	6.6	2.5	2	2	0	;
23	7.4	2.7	2	0	0.4056	;
24	0	3.1	2	0	0	;
25	0	2	2	0	0	;
26	6.8	5	2	0	0.4056]	;

nodes(:,3)=nodes(:,3)/max(nodes(:,2));
nodes(:,4)=nodes(:,4)/max(nodes(:,2));
nodes(:,2)=nodes(:,2)/max(nodes(:,2));





connectivity=[1	1	2	;
2	3	7	;
3	3	4	;
4	4	24	;
5	12	15	;
6	10	12	;
7	10	11	;
8	11	16	;
9	2	15	;
10	16	22	;
11	16	17	;
12	17	26	;
13	15	16	;
14	21	22	;
15	14	20	;
16	20	21	;
17	13	21	;
18	22	23	;
19	8	9	;
20	7	9	;
21	9	10	;
22	25	7	;
23	2	3	;
24	19	18	;
25	18	15	;
26	4	5	;
27	5	6	;
28	6	12	
]	;

figure
printNetwork(nodes,connectivity,20,'rete.pts',0,1);
