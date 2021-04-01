clear, clc, close all;

%%
% Приведение времен к одной длине (добавление нулей в интервале времен с 0 до 10 часов)
%%

% F = dir('*PS*.txt');
% for i = 1:length(F)
%     if ~F(i).isdir % если элемент НЕ папка
%     fname = F(i).name; % вытягиваем имя файла
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
%                 error('Неверное имя файла');
%             end
%         end
%     end
% end

clear, clc;
%%
% В имени файла дата начала эксперимента в формате "ГГММДД_ЧЧММСС"
%
% DataMIT - каждый столбец соответствует измерению PS волн
%           Строки:1    Время от начала эксперимента (с)
%                  2    Температура внутри камеры (°С)
%                  3    Поровое давление на входе (мВ)
%                  4    Поровое давление на выходе (мВ)
%                  5    Осевое давление (мВ)
%                  6    Боковое давление (мВ)
%                  7    [Пустой канал]
%                  8    Температура холодных концов термопары (°С)
%                  9    Размер образца (мм)
%                  10   Температура процесса (°С)
%                  11   Температура заданная (°С)
%                  12   Температура внутри термостата (°С)
% 
%   Time - вектор времени для осциллограмм (мкс)
%   DataP, DataS - каждый столбец - осциллограмма P/S волны (В)
% 
% Таким образом, чтобы получить все данные об измерении которое было
% произведено через 100 секунд от начала эксперимента необходимо в
% переменной DataMIT в первой строке найти ближайшее время, определить
% номер столбца и взять этот столбец из переменных DataMIT, DataP, DataS
%%

[~, Name, ~]= fileparts(cd);

F = dir('*PS_*.txt'); %Выбор данных по идентификатору

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
%         error('Время в процессе эксперимента было изменено!');
    end
    DataP(1:size(Data,1),count) = Data(1:size(Data,1),3);
    DataS(1:size(Data,1),count) = Data(1:size(Data,1),6);
    count = count+1;
end

% правка высоты по эталону
DataMIT(9,:) = DataMIT(9,:) - 0.458; % высота

save(strcat('..\',Name,'.mat'),'DataMIT','Time','DataP','DataS');
toc