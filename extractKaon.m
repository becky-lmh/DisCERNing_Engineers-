function kCharge = extractKaon(collision_tab)

% Remove particles that travel in negative z direction
collision_tab(collision_tab.p_z<0,:)=[];

% Extract pions
k_plus = collision_tab(ismember(collision_tab.Secondary,'kaon+'),:);
k_minus = collision_tab(ismember(collision_tab.Secondary,'kaon-'),:);

% Save pion data into a table
kCharge = [k_plus;k_minus];

% Remove pions with mass-energy < 200 MeV to ensure they will travel to all
% 3 detectors
kCharge(kCharge.mass_energy_MeV<200,:)=[];

% Remove pions that don't have pseudorapidity 1.38 < eta < 3.88
px = table2array(kCharge(:,"p_x"));
py = table2array(kCharge(:,"p_y"));
pz = table2array(kCharge(:,"p_z"));
eta = asinh(pz.*sin(atan(py./px))./py);
kCharge = [kCharge table(eta)];
kCharge((-1.38<kCharge.eta) & (kCharge.eta<1.38),:)=[];
kCharge(kCharge.eta>3.88,:)=[];
kCharge(kCharge.eta<-3.88,:)=[];



% k_plus = collision_tab(ismember(collision_tab.Secondary,'kaon+'),:);
% k_minus = collision_tab(ismember(collision_tab.Secondary,'kaon-'),:);
% proton = collision_tab(ismember(collision_tab.Secondary,'proton'),:);
% collision_charged = [pi_plus;pi_minus;k_plus;k_minus;proton];
% particles = categorical({'kaon+','kaon-','pi+','pi-','proton'});
% [C,ia,ic] = unique(table2array(collision_charged(:,2)));
% a_counts = accumarray(ic,1);
% bar(particles,a_counts);


end