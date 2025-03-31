%% ENTRYPOINT

disp('Please select which function to run:');
disp('1. (Centralized) Fully connectivity scenario');
disp('2. (Distributed) Robust MDS');
disp('3. (Centralized) MDS');


choice = input('Enter your choice (1, 2, 3 or 4): ');

switch choice
    case 1
        disp('Running Fully connectivity scenario...');
        main_emds_fully();
    case 2
        disp('Running Distributed scenario...');
        main_emds_distributed();
    case 3
        disp('Running Classical MDS...');
        main_classical_mds();    
    otherwise
        disp('Invalid choice. Please run the script again and select 1, 2, 3 or 4.');
end