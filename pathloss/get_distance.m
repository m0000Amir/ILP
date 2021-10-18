function d = get_distance(Ptr, Ltr, Gtr, Grecv, Lrecv, SOM, Precv, ...
    F, K)
% ������ ��������� ����� ����� ���������� � ������������
% 
% Ptr - �������� �����������,  [dBm]
% Ltr - ������ ������� � �������� ������, [dB]
% Gtr - ����������� �������� ���������� �������, [dBi]
% Lfs - ������ �������� � ��������� ������������, [dB]
% Grecv - ����������� �������� �������� �������,  [dBi]
% Lrecv - ������ ������� � �������� ������ � �������� ��������� �������, [dB]
% SOM - ����� �� ��������� ������� (System Operating Margin) [dB]
% Precv - ���������������� ���������

%% DEFAULT VALUE
% SOM = 10;
% Ltr = 0;
% Lrecv = 0;

%% LINK BUDGET
% Ptr - Ltr + Gtr - Lbf + Grecv - Lrecv = SOM + Precv

% free space path loss
Lfs = Ptr - Ltr + Gtr + Grecv - Lrecv - Precv - SOM;

%% PATH LOSS
% Lbf = 20lgF + 20lgD + K

% ��� �������, ���������� � ���, � ���������� � �����������, K = 92,45
% ��� �������, ���������� � ���, � ���������� � �����������, K = 32,4
% ��� �������, ���������� � ���, � ���������� � ������, K =  -27,55

% DISTANCE

d = floor(10^((Lfs - 20*log10(F) - K)/20));
end