% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function X = NTCIRSDmap(i,mapto)
%This functions maps the NTC pixels 1-128 into SRCCH format and vice versa
%See NTC Rosetta Stone: https://docs.google.com/spreadsheets/d/1RBM7bU40RPUYE5Yir00huaN4qAF6mQY3DwBZ1RnqrHs/edit#gid=1975841640

M=getmap();
switch mapto
    case {'pixel','pid'}
        %build 4D lookup table
        si = [12 4 4 8];
        spi = zeros(si);
        ind = sub2ind(si,M(:,1),M(:,2)+1,M(:,3)+1,M(:,4)+1);
        spi(ind) = 1:size(M,1); %populate locations with pixel numbers
        
        ind = sub2ind(si,i(:,1),i(:,2)+1,i(:,3)+1,i(:,4)+1); %find indices in 4D matrix
        X = spi(ind);
    case 'SRCCH'
        X = M(i,:);
end
end


function X=getmap()
%X=[...
% 12	0	1	6	0.00	1	1	8	2	2	12	6	15
% 12	0	1	4	0.00	1	1	7	2	2	12	5	13
% 12	0	1	2	0.00	1	1	6	2	2	12	5	11
% 12	0	1	0	0.00	1	1	5	2	2	12	5	9
% 12	0	2	6	0.00	1	1	4	3	3	12	5	23
% 12	0	2	4	0.00	1	1	3	3	3	12	5	21
% 12	0	2	2	0.00	1	1	2	3	3	12	5	19
% 12	0	2	0	0.00	1	1	1	3	3	12	5	17
% 12	0	1	7	0.00	1	2	8	2	2	12	5	16
% 12	0	1	5	0.00	1	2	7	2	2	12	5	14
% 12	0	1	3	0.00	1	2	6	2	2	12	5	12
% 12	0	1	1	0.00	1	2	5	2	2	12	5	10
% 12	0	2	7	0.00	1	2	4	3	3	12	5	24
% 12	0	2	5	0.00	1	2	3	3	3	12	5	22
% 12	0	2	3	0.00	1	2	2	3	3	12	5	20
% 12	0	2	1	0.00	1	2	1	3	3	12	5	18
% 12	0	0	1	0.00	1	3	8	1	1	12	5	2
% 12	0	0	3	0.00	1	3	7	1	1	12	5	4
% 12	0	0	5	0.00	1	3	6	1	1	12	5	6
% 12	0	0	7	0.00	1	3	5	1	1	12	5	8
% 12	0	3	1	0.00	1	3	4	4	4	12	5	26
% 12	0	3	3	0.00	1	3	3	4	4	12	5	28
% 12	0	3	5	0.00	1	3	2	4	4	12	5	30
% 12	0	3	7	0.00	1	3	1	4	4	12	5	32
% 12	0	0	0	0.00	1	4	8	1	1	12	5	1
% 12	0	0	2	0.00	1	4	7	1	1	12	5	3
% 12	0	0	4	0.00	1	4	6	1	1	12	5	5
% 12	0	0	6	0.00	1	4	5	1	1	12	5	7
% 12	0	3	0	0.00	1	4	4	4	4	12	5	25
% 12	0	3	2	0.00	1	4	3	4	4	12	5	27
% 12	0	3	4	0.00	1	4	2	4	4	12	5	29
% 12	0	3	6	0.00	1	4	1	4	4	12	5	31
% 12	1	3	6	0.00	1	5	8	8	8	12	5	63
% 12	1	3	4	0.00	1	5	7	8	8	12	5	61
% 12	1	3	2	0.00	1	5	6	8	8	12	5	59
% 12	1	3	0	0.00	1	5	5	8	8	12	5	57
% 12	1	0	6	0.00	1	5	4	5	5	12	5	39
% 12	1	0	4	0.00	1	5	3	5	5	12	5	37
% 12	1	0	2	0.00	1	5	2	5	5	12	5	35
% 12	1	0	0	0.00	1	5	1	5	5	12	5	33
% 12	1	3	7	0.00	1	6	8	8	8	12	5	64
% 12	1	3	5	0.00	1	6	7	8	8	12	5	62
% 12	1	3	3	0.00	1	6	6	8	8	12	5	60
% 12	1	3	1	0.00	1	6	5	8	8	12	5	58
% 12	1	0	7	0.00	1	6	4	5	5	12	5	40
% 12	1	0	5	0.00	1	6	3	5	5	12	5	38
% 12	1	0	3	0.00	1	6	2	5	5	12	5	36
% 12	1	0	1	0.00	1	6	1	5	5	12	5	34
% 12	1	2	1	0.00	1	7	8	7	7	12	5	50
% 12	1	2	3	0.00	1	7	7	7	7	12	5	52
% 12	1	2	5	0.00	1	7	6	7	7	12	5	54
% 12	1	2	7	0.00	1	7	5	7	7	12	5	56
% 12	1	1	1	0.00	1	7	4	6	6	12	5	42
% 12	1	1	3	0.00	1	7	3	6	6	12	5	44
% 12	1	1	5	0.00	1	7	2	6	6	12	5	46
% 12	1	1	7	0.00	1	7	1	6	6	12	5	48
% 12	1	2	0	0.00	1	8	8	7	7	12	5	49
% 12	1	2	2	0.00	1	8	7	7	7	12	5	51
% 12	1	2	4	0.00	1	8	6	7	7	12	5	53
% 12	1	2	6	0.00	1	8	5	7	7	12	5	55
% 12	1	1	0	0.00	1	8	4	6	6	12	5	41
% 12	1	1	2	0.00	1	8	3	6	6	12	5	43
% 12	1	1	4	0.00	1	8	2	6	6	12	5	45
% 12	1	1	6	0.00	1	8	1	6	6	12	5	47
% 12	3	1	6	0.00	2	8	1	14	14	12	6	111
% 12	3	1	4	0.00	2	8	2	14	14	12	6	109
% 12	3	1	2	0.00	2	8	3	14	14	12	6	107
% 12	3	1	0	0.00	2	8	4	14	14	12	6	105
% 12	3	2	6	0.00	2	8	5	15	15	12	6	119
% 12	3	2	4	0.00	2	8	6	15	15	12	6	117
% 12	3	2	2	0.00	2	8	7	15	15	12	6	115
% 12	3	2	0	0.00	2	8	8	15	15	12	6	113
% 12	3	1	7	0.00	2	7	1	14	14	12	6	112
% 12	3	1	5	0.00	2	7	2	14	14	12	6	110
% 12	3	1	3	0.00	2	7	3	14	14	12	6	108
% 12	3	1	1	0.00	2	7	4	14	14	12	6	106
% 12	3	2	7	0.00	2	7	5	15	15	12	6	120
% 12	3	2	5	0.00	2	7	6	15	15	12	6	118
% 12	3	2	3	0.00	2	7	7	15	15	12	6	116
% 12	3	2	1	0.00	2	7	8	15	15	12	6	114
% 12	3	0	1	0.00	2	6	1	13	13	12	6	98
% 12	3	0	3	0.00	2	6	2	13	13	12	6	100
% 12	3	0	5	0.00	2	6	3	13	13	12	6	102
% 12	3	0	7	0.00	2	6	4	13	13	12	6	104
% 12	3	3	1	0.00	2	6	5	16	16	12	6	122
% 12	3	3	3	0.00	2	6	6	16	16	12	6	124
% 12	3	3	5	0.00	2	6	7	16	16	12	6	126
% 12	3	3	7	0.00	2	6	8	16	16	12	6	128
% 12	3	0	0	0.00	2	5	1	13	13	12	6	97
% 12	3	0	2	0.00	2	5	2	13	13	12	6	99
% 12	3	0	4	0.00	2	5	3	13	13	12	6	101
% 12	3	0	6	0.00	2	5	4	13	13	12	6	103
% 12	3	3	0	0.00	2	5	5	16	16	12	6	121
% 12	3	3	2	0.00	2	5	6	16	16	12	6	123
% 12	3	3	4	0.00	2	5	7	16	16	12	6	125
% 12	3	3	6	0.00	2	5	8	16	16	12	6	127
% 12	2	3	6	0.00	2	4	1	12	12	12	6	95
% 12	2	3	4	0.00	2	4	2	12	12	12	6	93
% 12	2	3	2	0.00	2	4	3	12	12	12	6	91
% 12	2	3	0	0.00	2	4	4	12	12	12	6	89
% 12	2	0	6	0.00	2	4	5	9	9	12	6	71
% 12	2	0	4	0.00	2	4	6	9	9	12	6	69
% 12	2	0	2	0.00	2	4	7	9	9	12	6	67
% 12	2	0	0	0.00	2	4	8	9	9	12	6	65
% 12	2	3	7	0.00	2	3	1	12	12	12	6	96
% 12	2	3	5	0.00	2	3	2	12	12	12	6	94
% 12	2	3	3	0.00	2	3	3	12	12	12	6	92
% 12	2	3	1	0.00	2	3	4	12	12	12	6	90
% 12	2	0	7	0.00	2	3	5	9	9	12	6	72
% 12	2	0	5	0.00	2	3	6	9	9	12	6	70
% 12	2	0	3	0.00	2	3	7	9	9	12	6	68
% 12	2	0	1	0.00	2	3	8	9	9	12	6	66
% 12	2	2	1	0.00	2	2	1	11	11	12	6	82
% 12	2	2	3	0.00	2	2	2	11	11	12	6	84
% 12	2	2	5	0.00	2	2	3	11	11	12	6	86
% 12	2	2	7	0.00	2	2	4	11	11	12	6	88
% 12	2	1	1	0.00	2	2	5	10	10	12	6	74
% 12	2	1	3	0.00	2	2	6	10	10	12	6	76
% 12	2	1	5	0.00	2	2	7	10	10	12	6	78
% 12	2	1	7	0.00	2	2	8	10	10	12	6	80
% 12	2	2	0	0.00	2	1	1	11	11	12	6	81
% 12	2	2	2	0.00	2	1	2	11	11	12	6	83
% 12	2	2	4	0.00	2	1	3	11	11	12	6	85
% 12	2	2	6	0.00	2	1	4	11	11	12	6	87
% 12	2	1	0	0.00	2	1	5	10	10	12	6	73
% 12	2	1	2	0.00	2	1	6	10	10	12	6	75
% 12	2	1	4	0.00	2	1	7	10	10	12	6	77
% 12	2	1	6	0.00	2	1	8	10	10	12	6	79]; %Omron2 Mapping
% 
% return;

% %NEW ArrayX Card Mapping from 7/2017
% X=[...
% 12	0	0	0	0.00	1	1	8	2	2	12	6	15
% 12	0	0	1	0.00	1	1	7	2	2	12	5	13
% 12	0	0	2	0.00	1	1	6	2	2	12	5	11
% 12	0	0	3	0.00	1	1	5	2	2	12	5	9
% 12	0	1	0	0.00	1	1	4	3	3	12	5	23
% 12	0	1	1	0.00	1	1	3	3	3	12	5	21
% 12	0	1	3	0.00	1	1	2	3	3	12	5	19
% 12	0	1	2	0.00	1	1	1	3	3	12	5	17
% 12	0	0	6	0.00	1	2	8	2	2	12	5	16
% 12	0	0	7	0.00	1	2	7	2	2	12	5	14
% 12	0	0	5	0.00	1	2	6	2	2	12	5	12
% 12	0	0	4	0.00	1	2	5	2	2	12	5	10
% 12	0	1	7	0.00	1	2	4	3	3	12	5	24
% 12	0	1	6	0.00	1	2	3	3	3	12	5	22
% 12	0	1	5	0.00	1	2	2	3	3	12	5	20
% 12	0	1	4	0.00	1	2	1	3	3	12	5	18
% 12	1	0	2	0.00	1	3	8	1	1	12	5	2
% 12	1	0	0	0.00	1	3	7	1	1	12	5	4
% 12	1	0	7	0.00	1	3	6	1	1	12	5	6
% 12	1	0	3	0.00	1	3	5	1	1	12	5	8
% 12	1	1	3	0.00	1	3	4	4	4	12	5	26
% 12	1	1	4	0.00	1	3	3	4	4	12	5	28
% 12	1	1	2	0.00	1	3	2	4	4	12	5	30
% 12	1	1	0	0.00	1	3	1	4	4	12	5	32
% 12	1	0	4	0.00	1	4	8	1	1	12	5	1
% 12	1	0	5	0.00	1	4	7	1	1	12	5	3
% 12	1	0	6	0.00	1	4	6	1	1	12	5	5
% 12	1	0	1	0.00	1	4	5	1	1	12	5	7
% 12	1	1	7	0.00	1	4	4	4	4	12	5	25
% 12	1	1	5	0.00	1	4	3	4	4	12	5	27
% 12	1	1	1	0.00	1	4	2	4	4	12	5	29
% 12	1	1	6	0.00	1	4	1	4	4	12	5	31
% 12	2	0	0	0.00	1	5	8	8	8	12	5	63
% 12	2	0	4	0.00	1	5	7	8	8	12	5	61
% 12	2	0	3	0.00	1	5	6	8	8	12	5	59
% 12	2	1	0	0.00	1	5	5	8	8	12	5	57
% 12	2	1	7	0.00	1	5	4	5	5	12	5	39
% 12	2	1	6	0.00	1	5	3	5	5	12	5	37
% 12	2	1	5	0.00	1	5	2	5	5	12	5	35
% 12	2	1	3	0.00	1	5	1	5	5	12	5	33
% 12	2	0	7	0.00	1	6	8	8	8	12	5	64
% 12	2	0	6	0.00	1	6	7	8	8	12	5	62
% 12	2	0	1	0.00	1	6	6	8	8	12	5	60
% 12	2	0	2	0.00	1	6	5	8	8	12	5	58
% 12	2	0	5	0.00	1	6	4	5	5	12	5	40
% 12	2	1	1	0.00	1	6	3	5	5	12	5	38
% 12	2	1	2	0.00	1	6	2	5	5	12	5	36
% 12	2	1	4	0.00	1	6	1	5	5	12	5	34
% 12	3	0	4	0.00	1	7	8	7	7	12	5	50
% 12	3	0	1	0.00	1	7	7	7	7	12	5	52
% 12	3	0	0	0.00	1	7	6	7	7	12	5	54
% 12	3	0	5	0.00	1	7	5	7	7	12	5	56
% 12	3	1	3	0.00	1	7	4	6	6	12	5	42
% 12	3	1	2	0.00	1	7	3	6	6	12	5	44
% 12	3	1	6	0.00	1	7	2	6	6	12	5	46
% 12	3	1	7	0.00	1	7	1	6	6	12	5	48
% 12	3	0	3	0.00	1	8	8	7	7	12	5	49
% 12	3	0	2	0.00	1	8	7	7	7	12	5	51
% 12	3	0	6	0.00	1	8	6	7	7	12	5	53
% 12	3	0	7	0.00	1	8	5	7	7	12	5	55
% 12	3	1	1	0.00	1	8	4	6	6	12	5	41
% 12	3	1	4	0.00	1	8	3	6	6	12	5	43
% 12	3	1	5	0.00	1	8	2	6	6	12	5	45
% 12	3	1	0	0.00	1	8	1	6	6	12	5	47
% 12	3	3	0	0.00	2	8	1	14	14	12	6	111
% 12	3	3	5	0.00	2	8	2	14	14	12	6	109
% 12	3	3	4	0.00	2	8	3	14	14	12	6	107
% 12	3	3	1	0.00	2	8	4	14	14	12	6	105
% 12	3	2	7	0.00	2	8	5	15	15	12	6	119
% 12	3	2	6	0.00	2	8	6	15	15	12	6	117
% 12	3	2	2	0.00	2	8	7	15	15	12	6	115
% 12	3	2	3	0.00	2	8	8	15	15	12	6	113
% 12	3	3	7	0.00	2	7	1	14	14	12	6	112
% 12	3	3	6	0.00	2	7	2	14	14	12	6	110
% 12	3	3	2	0.00	2	7	3	14	14	12	6	108
% 12	3	3	3	0.00	2	7	4	14	14	12	6	106
% 12	3	2	5	0.00	2	7	5	15	15	12	6	120
% 12	3	2	0	0.00	2	7	6	15	15	12	6	118
% 12	3	2	1	0.00	2	7	7	15	15	12	6	116
% 12	3	2	4	0.00	2	7	8	15	15	12	6	114
% 12	2	3	4	0.00	2	6	1	13	13	12	6	98
% 12	2	3	2	0.00	2	6	2	13	13	12	6	100
% 12	2	3	1	0.00	2	6	3	13	13	12	6	102
% 12	2	2	5	0.00	2	6	4	13	13	12	6	104
% 12	2	2	2	0.00	2	6	5	16	16	12	6	122
% 12	2	2	1	0.00	2	6	6	16	16	12	6	124
% 12	2	2	6	0.00	2	6	7	16	16	12	6	126
% 12	2	2	7	0.00	2	6	8	16	16	12	6	128
% 12	2	3	3	0.00	2	5	1	13	13	12	6	97
% 12	2	3	5	0.00	2	5	2	13	13	12	6	99
% 12	2	3	6	0.00	2	5	3	13	13	12	6	101
% 12	2	3	7	0.00	2	5	4	13	13	12	6	103
% 12	2	3	0	0.00	2	5	5	16	16	12	6	121
% 12	2	2	3	0.00	2	5	6	16	16	12	6	123
% 12	2	2	4	0.00	2	5	7	16	16	12	6	125
% 12	2	2	0	0.00	2	5	8	16	16	12	6	127
% 12	1	3	6	0.00	2	4	1	12	12	12	6	95
% 12	1	3	1	0.00	2	4	2	12	12	12	6	93
% 12	1	3	5	0.00	2	4	3	12	12	12	6	91
% 12	1	3	7	0.00	2	4	4	12	12	12	6	89
% 12	1	2	1	0.00	2	4	5	9	9	12	6	71
% 12	1	2	6	0.00	2	4	6	9	9	12	6	69
% 12	1	2	5	0.00	2	4	7	9	9	12	6	67
% 12	1	2	4	0.00	2	4	8	9	9	12	6	65
% 12	1	3	0	0.00	2	3	1	12	12	12	6	96
% 12	1	3	2	0.00	2	3	2	12	12	12	6	94
% 12	1	3	4	0.00	2	3	3	12	12	12	6	92
% 12	1	3	3	0.00	2	3	4	12	12	12	6	90
% 12	1	2	3	0.00	2	3	5	9	9	12	6	72
% 12	1	2	7	0.00	2	3	6	9	9	12	6	70
% 12	1	2	0	0.00	2	3	7	9	9	12	6	68
% 12	1	2	2	0.00	2	3	8	9	9	12	6	66
% 12	0	3	4	0.00	2	2	1	11	11	12	6	82
% 12	0	3	5	0.00	2	2	2	11	11	12	6	84
% 12	0	3	6	0.00	2	2	3	11	11	12	6	86
% 12	0	3	7	0.00	2	2	4	11	11	12	6	88
% 12	0	2	4	0.00	2	2	5	10	10	12	6	74
% 12	0	2	5	0.00	2	2	6	10	10	12	6	76
% 12	0	2	7	0.00	2	2	7	10	10	12	6	78
% 12	0	2	6	0.00	2	2	8	10	10	12	6	80
% 12	0	3	2	0.00	2	1	1	11	11	12	6	81
% 12	0	3	3	0.00	2	1	2	11	11	12	6	83
% 12	0	3	1	0.00	2	1	3	11	11	12	6	85
% 12	0	3	0	0.00	2	1	4	11	11	12	6	87
% 12	0	2	3	0.00	2	1	5	10	10	12	6	73
% 12	0	2	2	0.00	2	1	6	10	10	12	6	75
% 12	0	2	1	0.00	2	1	7	10	10	12	6	77
% 12	0	2	0	0.00	2	1	8	10	10	12	6	79];

%UPDATED ArrayX Card Mapping from 10/13/2017
X=[...
12	0	0	0	0.00	1	nan	nan	nan	nan	nan	nan	15
12	0	0	1	0.00	1	nan	nan	nan	nan	nan	nan	13
12	0	0	2	0.00	1	nan	nan	nan	nan	nan	nan	11
12	0	0	3	0.00	1	nan	nan	nan	nan	nan	nan	9
12	0	1	0	0.00	1	nan	nan	nan	nan	nan	nan	23
12	0	1	1	0.00	1	nan	nan	nan	nan	nan	nan	21
12	0	1	3	0.00	1	nan	nan	nan	nan	nan	nan	19
12	0	1	2	0.00	1	nan	nan	nan	nan	nan	nan	17
12	0	0	6	0.00	1	nan	nan	nan	nan	nan	nan	16
12	0	0	7	0.00	1	nan	nan	nan	nan	nan	nan	14
12	0	0	5	0.00	1	nan	nan	nan	nan	nan	nan	12
12	0	0	4	0.00	1	nan	nan	nan	nan	nan	nan	10
12	0	1	7	0.00	1	nan	nan	nan	nan	nan	nan	24
12	0	1	6	0.00	1	nan	nan	nan	nan	nan	nan	22
12	0	1	5	0.00	1	nan	nan	nan	nan	nan	nan	20
12	0	1	4	0.00	1	nan	nan	nan	nan	nan	nan	18
12	1	0	2	0.00	1	nan	nan	nan	nan	nan	nan	2
12	1	0	0	0.00	1	nan	nan	nan	nan	nan	nan	4
12	1	0	7	0.00	1	nan	nan	nan	nan	nan	nan	6
12	1	0	3	0.00	1	nan	nan	nan	nan	nan	nan	8
12	1	1	3	0.00	1	nan	nan	nan	nan	nan	nan	26
12	1	1	4	0.00	1	nan	nan	nan	nan	nan	nan	28
12	1	1	2	0.00	1	nan	nan	nan	nan	nan	nan	30
12	1	1	0	0.00	1	nan	nan	nan	nan	nan	nan	32
12	1	0	4	0.00	1	nan	nan	nan	nan	nan	nan	1
12	1	0	5	0.00	1	nan	nan	nan	nan	nan	nan	3
12	1	0	6	0.00	1	nan	nan	nan	nan	nan	nan	5
12	1	0	1	0.00	1	nan	nan	nan	nan	nan	nan	7
12	1	1	7	0.00	1	nan	nan	nan	nan	nan	nan	25
12	1	1	5	0.00	1	nan	nan	nan	nan	nan	nan	27
12	1	1	1	0.00	1	nan	nan	nan	nan	nan	nan	29
12	1	1	6	0.00	1	nan	nan	nan	nan	nan	nan	31
12	2	0	0	0.00	1	nan	nan	nan	nan	nan	nan	63
12	2	0	4	0.00	1	nan	nan	nan	nan	nan	nan	61
12	2	0	3	0.00	1	nan	nan	nan	nan	nan	nan	59
12	2	1	0	0.00	1	nan	nan	nan	nan	nan	nan	57
12	2	1	7	0.00	1	nan	nan	nan	nan	nan	nan	39
12	2	1	6	0.00	1	nan	nan	nan	nan	nan	nan	37
12	2	1	5	0.00	1	nan	nan	nan	nan	nan	nan	35
12	2	1	3	0.00	1	nan	nan	nan	nan	nan	nan	33
12	2	0	7	0.00	1	nan	nan	nan	nan	nan	nan	64
12	2	0	6	0.00	1	nan	nan	nan	nan	nan	nan	62
12	2	0	1	0.00	1	nan	nan	nan	nan	nan	nan	60
12	2	0	2	0.00	1	nan	nan	nan	nan	nan	nan	58
12	2	0	5	0.00	1	nan	nan	nan	nan	nan	nan	40
12	2	1	1	0.00	1	nan	nan	nan	nan	nan	nan	38
12	2	1	2	0.00	1	nan	nan	nan	nan	nan	nan	36
12	2	1	4	0.00	1	nan	nan	nan	nan	nan	nan	34
12	3	0	4	0.00	1	nan	nan	nan	nan	nan	nan	50
12	3	0	1	0.00	1	nan	nan	nan	nan	nan	nan	52
12	3	0	0	0.00	1	nan	nan	nan	nan	nan	nan	54
12	3	0	5	0.00	1	nan	nan	nan	nan	nan	nan	56
12	3	1	3	0.00	1	nan	nan	nan	nan	nan	nan	42
12	3	1	2	0.00	1	nan	nan	nan	nan	nan	nan	44
12	3	1	6	0.00	1	nan	nan	nan	nan	nan	nan	46
12	3	1	7	0.00	1	nan	nan	nan	nan	nan	nan	48
12	3	0	3	0.00	1	nan	nan	nan	nan	nan	nan	49
12	3	0	2	0.00	1	nan	nan	nan	nan	nan	nan	51
12	3	0	6	0.00	1	nan	nan	nan	nan	nan	nan	53
12	3	0	7	0.00	1	nan	nan	nan	nan	nan	nan	55
12	3	1	1	0.00	1	nan	nan	nan	nan	nan	nan	41
12	3	1	4	0.00	1	nan	nan	nan	nan	nan	nan	43
12	3	1	5	0.00	1	nan	nan	nan	nan	nan	nan	45
12	3	1	0	0.00	1	nan	nan	nan	nan	nan	nan	47
12	3	3	2	0.00	2	nan	nan	nan	nan	nan	nan	111
12	3	3	0	0.00	2	nan	nan	nan	nan	nan	nan	109
12	3	3	4	0.00	2	nan	nan	nan	nan	nan	nan	107
12	3	3	6	0.00	2	nan	nan	nan	nan	nan	nan	105
12	3	2	0	0.00	2	nan	nan	nan	nan	nan	nan	119
12	3	2	2	0.00	2	nan	nan	nan	nan	nan	nan	117
12	3	2	4	0.00	2	nan	nan	nan	nan	nan	nan	115
12	3	2	6	0.00	2	nan	nan	nan	nan	nan	nan	113
12	3	3	1	0.00	2	nan	nan	nan	nan	nan	nan	112
12	3	3	3	0.00	2	nan	nan	nan	nan	nan	nan	110
12	3	3	5	0.00	2	nan	nan	nan	nan	nan	nan	108
12	3	3	7	0.00	2	nan	nan	nan	nan	nan	nan	106
12	3	2	1	0.00	2	nan	nan	nan	nan	nan	nan	120
12	3	2	3	0.00	2	nan	nan	nan	nan	nan	nan	118
12	3	2	7	0.00	2	nan	nan	nan	nan	nan	nan	116
12	3	2	5	0.00	2	nan	nan	nan	nan	nan	nan	114
12	2	3	6	0.00	2	nan	nan	nan	nan	nan	nan	98
12	2	3	2	0.00	2	nan	nan	nan	nan	nan	nan	100
12	2	3	1	0.00	2	nan	nan	nan	nan	nan	nan	102
12	2	3	0	0.00	2	nan	nan	nan	nan	nan	nan	104
12	2	2	0	0.00	2	nan	nan	nan	nan	nan	nan	122
12	2	2	7	0.00	2	nan	nan	nan	nan	nan	nan	124
12	2	2	6	0.00	2	nan	nan	nan	nan	nan	nan	126
12	2	2	2	0.00	2	nan	nan	nan	nan	nan	nan	128
12	2	3	5	0.00	2	nan	nan	nan	nan	nan	nan	97
12	2	3	4	0.00	2	nan	nan	nan	nan	nan	nan	99
12	2	3	3	0.00	2	nan	nan	nan	nan	nan	nan	101
12	2	3	7	0.00	2	nan	nan	nan	nan	nan	nan	103
12	2	2	4	0.00	2	nan	nan	nan	nan	nan	nan	121
12	2	2	5	0.00	2	nan	nan	nan	nan	nan	nan	123
12	2	2	3	0.00	2	nan	nan	nan	nan	nan	nan	125
12	2	2	1	0.00	2	nan	nan	nan	nan	nan	nan	127
12	1	3	0	0.00	2	nan	nan	nan	nan	nan	nan	95
12	1	3	3	1.00	2	nan	nan	nan	nan	nan	nan	93
12	1	3	5	1.00	2	nan	nan	nan	nan	nan	nan	91
12	1	3	7	0.00	2	nan	nan	nan	nan	nan	nan	89
12	1	3	6	0.00	2	nan	nan	nan	nan	nan	nan	71
12	1	2	0	0.00	2	nan	nan	nan	nan	nan	nan	69
12	1	2	1	0.00	2	nan	nan	nan	nan	nan	nan	67
12	1	2	6	0.00	2	nan	nan	nan	nan	nan	nan	65
12	1	3	1	0.00	2	nan	nan	nan	nan	nan	nan	96
12	1	3	2	1.00	2	nan	nan	nan	nan	nan	nan	94
12	1	3	4	0.00	2	nan	nan	nan	nan	nan	nan	92
12	1	2	3	1.00	2	nan	nan	nan	nan	nan	nan	90
12	1	2	2	0.00	2	nan	nan	nan	nan	nan	nan	72
12	1	2	4	0.00	2	nan	nan	nan	nan	nan	nan	70
12	1	2	5	0.00	2	nan	nan	nan	nan	nan	nan	68
12	1	2	7	0.00	2	nan	nan	nan	nan	nan	nan	66
12	0	3	7	0.00	2	nan	nan	nan	nan	nan	nan	82
12	0	3	5	0.00	2	nan	nan	nan	nan	nan	nan	84
12	0	3	2	0.00	2	nan	nan	nan	nan	nan	nan	86
12	0	3	0	0.00	2	nan	nan	nan	nan	nan	nan	88
12	0	2	3	0.00	2	nan	nan	nan	nan	nan	nan	74
12	0	2	6	0.00	2	nan	nan	nan	nan	nan	nan	76
12	0	2	4	0.00	2	nan	nan	nan	nan	nan	nan	78
12	0	2	1	0.00	2	nan	nan	nan	nan	nan	nan	80
12	0	3	6	0.00	2	nan	nan	nan	nan	nan	nan	81
12	0	3	3	0.00	2	nan	nan	nan	nan	nan	nan	83
12	0	3	1	0.00	2	nan	nan	nan	nan	nan	nan	85
12	0	3	4	0.00	2	nan	nan	nan	nan	nan	nan	87
12	0	2	7	0.00	2	nan	nan	nan	nan	nan	nan	73
12	0	2	5	0.00	2	nan	nan	nan	nan	nan	nan	75
12	0	2	2	0.00	2	nan	nan	nan	nan	nan	nan	77
12	0	2	0	0.00	2	nan	nan	nan	nan	nan	nan	79];
end