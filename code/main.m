%% ENTRYPOINT

disp('Please select which function to run:');
disp('1. (Centralized) Partial connectivity scenario');
disp('2. (Centralized) Fully connectivity scenario');
disp('3. (Distributed) Robust MDS');
disp('4. (Centralized) MDS');


choice = input('Enter your choice (1, 2, 3 or 4): ');

switch choice
    case 1
        disp('Running Partial connectivity scenario...');
        main_emds_partial();
    case 2
        disp('Running Fully connectivity scenario...');
        main_emds_fully();
    case 3
        disp('Running Distributed scenario...');
        main_emds_distributed();
    case 4
        disp('Running Classical MDS...');
        main_classical_mds();
    
    otherwise
        disp('Invalid choice. Please run the script again and select 1, 2, 3 or 4.');
end