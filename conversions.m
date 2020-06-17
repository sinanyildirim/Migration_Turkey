Iller = {'Adana-1',
'Adiyaman-2',
'Afyonkarahisar-3'
'Agri-4',
'Aksaray-68',
'Amasya-5',
'Ankara-6',
'Antalya-7',
'Ardahan-75',
'Artvin-8',
'Aydin-9',
'Balikesir-10',
'Bartin-74',
'Batman-72',
'Bayburt-69',
'Bilecik-11',
'Bingol-12',
'Bitlis-13',
'Bolu-14',
'Burdur-15',
'Bursa-16',
'Canakkale-17',
'Cankiri-18',
'Corum-19',
'Denizli-20',
'Diyarbakir-21',
'Duzce-81',
'Edirne-22',
'Elazig-23',
'Erzincan-24',
'Erzurum-25',
'Eskisehir-26',
'Gaziantep-27',
'Giresun-28',
'Gumushane-29',
'Hakkari-30',
'Hatay-31',
'Igdir-76',
'Isparta-32',
'Istanbul-34',
'Izmir-35',
'Kahramanmaras-46',
'Karabuk-78',
'Karaman-70',
'Kars-36',
'Kastamonu-37',
'Kayseri-38',
'Kilis-79',
'Kirikkale-71',
'Kirklareli-39',
'Kirsehir-40',
'Kocaeli-41',
'Konya-42',
'Kutahya-43',
'Malatya-44',
'Manisa-45',
'Mardin-47',
'Mersin-33',
'Mugla-48',
'Mus-49',
'Nevsehir-50',
'Nigde-51',
'Ordu-52',
'Osmaniye-80',
'Rize-53',
'Sakarya-54',
'Samsun-55',
'Sanliurfa-63',
'Siirt-56',
'Sinop-57',
'Sirnak-73',
'Sivas-58',
'Tekirdag-59',
'Tokat-60',
'Trabzon-61',
'Tunceli-62',
'Usak-64',
'Van-65',
'Yalova-77',
'Yozgat-66',
'Zonguldak-67'};

load('ilmeseafe_data.mat');

Duzey_isim = {'TR10'
    'TR21'
    'TR22'
    'TR31'
    'TR32'
    'TR33'
    'TR41'
    'TR42'
    'TR51'
    'TR52'
    'TR61'
    'TR62'
    'TR63'
    'TR71'
    'TR72'
    'TR81'
    'TR82'
    'TR83'
    'TR90'
    'TRA1'
    'TRA2'
    'TRB1'
    'TRB2'
    'TRC1'
    'TRC2'
    'TRC3'};

Duzey_2_Il_abc = {[40], [73, 28, 50], [12, 22], [41], [11, 25, 59], [56, 3, 54, 77], ...
    [16, 21, 32], [52, 66, 27, 19, 79], [7], [53, 44], [8, 39, 20], ...
    [1, 58], [37, 42, 64], [49, 5, 62, 61, 51], [47, 72, 80],....
    [81, 43, 13], [46, 23, 70], [67, 74, 24, 6], [75, 63, 34, 65, 10, 35],...
    [31, 30, 15], [4, 45, 38, 9], [55, 29, 17, 76], [78, 60, 18, 36],...
    [33, 2, 48], [68, 26], [57, 14, 69, 71]};

% from cities (in alphabetical order) to nuts-2
Il_abc_2_Duzey = zeros(81, 1);
for i = 1:81
    for j = 1:26;
        if find(Duzey_2_Il_abc{j} == i)
            Il_abc_2_Duzey(i) = j;
            break;
        end
    end
end            

% from plaka (in numerical order) to nuts-2
il_plaka_2_abc = zeros(81, 1);
for i = 1:81
    for j = 1:81
        if i == il_abc_2_plaka(j)
            il_plaka_2_abc(i) = j;
            break;
        end
    end
end



Il_plaka_2_Duzey = Il_abc_2_Duzey(il_plaka_2_abc);

load('bolge_adlari.mat');
load('bolge_adlari2.mat');

bolge_adlari_1_to_2 = zeros(1, 26);
for i = 1:26
    for j = 1:26
        if strcmp(bolge_adlari2{j}, bolge_adlari{i}) == 1
            bolge_adlari_1_to_2(i) = j;
            break;
        end
    end
end

clear i, j;

save('conversions');

% 
