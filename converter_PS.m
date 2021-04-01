clear, clc, close all;

%%
% ���������� ������ � ����� ����� (���������� ����� � ��������� ������ � 0 �� 10 �����)
%%

% F = dir('*PS*.txt');
% for i = 1:length(F)
%     if ~F(i).isdir % ���� ������� �� �����
%     fname = F(i).name; % ���������� ��� �����
%         if fname(16)=='_'
%             new_fname = fname;
%             str = sprintf('%d| File: %s -> %s',i ,fname, new_fname);
%             disp(str);
%         else
%             if fname(15)=='_'
%                 new_fname = strcat(fname(1:9),'0',fname(10:end));
%                 str = sprintf('%d| File: %s -> %s',i ,fname, new_fname);
%                 disp(str);
%                 movefile(F(i,1).name,new_fname);        %# Rename the file
%             else
%                 error('�������� ��� �����');
%             end
%         end
%     end
% end

clear, clc;
%%
% � ����� ����� ���� ������ ������������ � ������� "������_������"
%
% DataMIT - ������ ������� ������������� ��������� PS ����
%           ������:1    ����� �� ������ ������������ (�)
%                  2    ����������� ������ ������ (��)
%                  3    ������� �������� �� ����� (��)
%                  4    ������� �������� �� ������ (��)
%                  5    ������ �������� (��)
%                  6    ������� �������� (��)
%                  7    [������ �����]
%                  8    ����������� �������� ������ ��������� (��)
%                  9    ������ ������� (��)
%                  10   ����������� �������� (��)
%                  11   ����������� �������� (��)
%                  12   ����������� ������ ���������� (��)
% 
%   Time - ������ ������� ��� ������������ (���)
%   DataP, DataS - ������ ������� - ������������� P/S ����� (�)
% 
% ����� �������, ����� �������� ��� ������ �� ��������� ������� ����
% ����������� ����� 100 ������ �� ������ ������������ ���������� �
% ���������� DataMIT � ������ ������ ����� ��������� �����, ����������
% ����� ������� � ����� ���� ������� �� ���������� DataMIT, DataP, DataS
%%

[~, Name, ~]= fileparts(cd);

F = dir('*PS_*.txt'); %����� ������ �� ��������������

step = 2;
copyFlag = true;

tic
count = 1; 
fileID = 1:step:length(F);
DataMIT = zeros(16,numel(fileID));
for i = fileID
    fprintf('Trace ''%d'' from ''%d''.\n', i, length(F));
    fname = F(i).name;
    
%     fid = fopen(fname);
%     Q = textscan(fid,'%5s',1,'delimiter','\r')
%     fclose(fid);
    
    Data = importdata(fname,'\t',0);
    DataMIT(:,count) = (Data(1,:))';
    
    DataT = importdata(fname,'\t',15);
    Data = DataT.data;
    Data(1,:) = [];
    
    k = 1;
    while  Data(k+1,1) < 0
        k = k + 1;
    end
    Data(1:k,:) = [];
    
    if (i == 1)
        Time = Data(:,1);
        DataP = zeros(size(Data,1),numel(fileID));
        DataS = zeros(size(Data,1),numel(fileID));
    end
    
    if (Time(1:size(Data,1),1) == Data(1:size(Data,1),1)*1000000)
        Time(1:size(Data,1),1) = Data(1:size(Data,1),1)*1000000;
    else
%         error('����� � �������� ������������ ���� ��������!');
    end
    DataP(1:size(Data,1),count) = Data(1:size(Data,1),3);
    DataS(1:size(Data,1),count) = Data(1:size(Data,1),6);
    count = count+1;
end

% ������ ������ �� �������
DataMIT(9,:) = DataMIT(9,:) - 0.458; % ������

save(strcat('..\',Name,'.mat'),'DataMIT','Time','DataP','DataS');
toc