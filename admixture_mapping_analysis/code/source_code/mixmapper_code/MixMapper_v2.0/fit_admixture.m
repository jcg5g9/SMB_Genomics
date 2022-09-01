% FIT_ADMIXTURE Helper function that creates a scaffold tree using
% neighbor joining and finds the best-fit placement of one or two
% additional populations as two- or three-way admixtures.
% 
%   [TREE,FIT] = FIT_ADMIXTURE(POP_DATA,SCAFFOLD_POP_NAMES,ADMIXED1)
%   computes best-fit parameters for fitting population ADMIXED1 as a
%   mixture of two populations branching from the scaffold tree.  POP_DATA
%   is a structure containing the following fields:
%   - f2_all (n x n matrix of pairwise f2 statistics)
%   - het (n vector of heterozygosity values)
%   - pop_names (n cell array)
%   SCAFFOLD_POP_NAMES is a cell array containing the populations to use
%   in the scaffold.
%   TREE is a phytree object (from the Matlab Bioinformatics Toolbox).
%   FIT is a structure containing parameters of the fit.
%
%   [TREE,FIT] = FIT_ADMIXTURE(...,ADMIXED1,ADMIXED2) fits ADMIXED1 as a
%   two-way admixture as above and subsequently fits ADMIXED2 as a mixture
%   of the ADMIXED1 branch and a third branch of the scaffold tree.
%
%   [TREE,FIT] = FIT_ADMIXTURE(...,ADMIXED1,ADMIXED2,OPTIONS) allows
%   setting additional options by passing a structure OPTIONS with the
%   following possible fields:
%   - drift_units: set to 1 to return branch lengths on the tree in drift
%     units (i.e., normalized by heterozygosity in the parent); note that as
%     drift units are not additive, the tree is made by first fitting in f2
%     units and then converting each branch length
%   - print_output: set to 0 to suppress output
%   - no_opt3way: set to 1 to revert to pre-2.0 3-way mixture fitting
%   - test_2v3: eliminate last equation from 3-way fitting for fair
%     comparison between 2-way and 3-way models
%
%   TREE = FIT_ADMIXTURE(POP_DATA,SCAFFOLD_POP_NAMES) creates and returns
%   only the scaffold tree.


% Mark Lipson and Po-Ru Loh, 5/29/2014

function [tree,fit] = fit_admixture(pop_data,scaffold_pop_names,admixed1,admixed2,options)

f2_all = pop_data.f2;
het = pop_data.h;
pop_names = pop_data.pop_names;

if numel(unique(pop_names)) ~= numel(pop_names)
    error('pop_data.pop_names contains duplicate(s)')
end
if numel(unique(scaffold_pop_names)) ~= numel(scaffold_pop_names)
    error('scaffold_pop_names contains duplicate(s)')
end

if nargin < 3
    admixed1 = '';
end
check_admixed_pop_valid(admixed1)

if nargin < 4
    admixed2 = '';
end
check_admixed_pop_valid(admixed2)

if nargin < 5
    options = [];
end

if check_option(options,'print_output')
    fid = 1; % print to screen
else % don't print output
    fid = fopen('/dev/null','w');
    if fid == -1
        fid = fopen('ignore_output.txt','w');
        %fid = 1; % give up and just display output on screen
    end
end

COEFFS_MULT = 1000; % multiplier on coeffs to rescale resnorms close to 1 for optimization

%%%%%%%%%% make tree from scaffold populations %%%%%%%%%%%

num_scaffold = length(scaffold_pop_names);
num_edges = 2*(num_scaffold-1);
num_nodes = num_edges+1;
is_scaffold = false(length(pop_names),1);
for p = 1:num_scaffold
    is_scaffold = is_scaffold | strcmp(pop_names,scaffold_pop_names{p});
    if sum(is_scaffold) ~= p
        error('population scaffold_pop_names{%d} = %s not found in pop_data.pop_names; check spelling?',p,scaffold_pop_names{p})
    end
end
tree = seqneighjoin(f2_all(is_scaffold,is_scaffold),'equivar',pop_names(is_scaffold));
node_names = get(tree,'nodenames');
% beware: scaffold_pop_names and node_names{1:num_scaffold} do not have same order!

%%%%%%%%%% find best root location according to heterozygosity data %%%%%%%%%%%

best_resnorm = 1e9;
best_branch = 1;
edge_lens = get(tree,'distances');
for branch = 1:num_nodes-1 % not root
    % reroot at distance very close to branch node to disambiguate which
    % child node corresponds to the original branch node in the old tree
    % (for use in rerooting with the real distance later)
    [resnorm,~,d] = ...
        fit_het(f2_all,het,pop_names,scaffold_pop_names,reroot(tree,getbyname(tree,node_names{branch},'exact',true),1e-9));
    if resnorm < best_resnorm && d < edge_lens(branch)
        best_resnorm = resnorm;
        best_branch = branch;
        reroot_dist = d;
    end
end

if best_resnorm < 1e9
    tree = reroot(tree,getbyname(tree,node_names{best_branch},'exact',true),reroot_dist);
end
% topology may have changed from previous fit_het upon reroot with different distance
[resnorm,het_nodes] = fit_het(f2_all,het,pop_names,scaffold_pop_names,tree);
fit.het_nodes = het_nodes;

tree_dist = pdist(tree,'nodes','all','squareform',true);

%%%%%%%%%% create table with parent of each node %%%%%%%%%%%

pointers = get(tree,'pointers');
parents = zeros(length(node_names),1);
for ptr = 1:size(pointers,1)
    for child=1:2
        parents(pointers(ptr,child)) = ptr + num_scaffold;
    end
end

%%%%%%%%%% create tree to return with distances in drift units %%%%%%%%%%%

node_names = get(tree,'nodenames');
edge_lens = get(tree,'distances');
drift_lens = edge_lens;
if check_option(options,'drift_units')
    for edge = 1:num_edges
        drift_lens(edge) = 2*edge_lens(edge) / het_nodes(parents(edge));
        %        fprintf('Branch %s: %f = 2*%f / %f\n',node_names{edge},drift_lens(edge),edge_lens(edge),het_nodes(parents(edge)));
    end
end
tree = phytree(pointers,drift_lens,node_names);
edge_lens = edge_lens(1:num_edges); % get rid of trailing 0 (root)

if isempty(admixed1)
    return
end

%%%%%%%%%% create lookup of scaffold_pop indices -> indices into f2 array %%%%%%%%%%%

scaffold_f2_inds = zeros(num_scaffold,1);
for p = 1:num_scaffold
    scaffold_f2_inds(p) = find_str(node_names{p});
end
admixed_ind1 = find_str(admixed1);
admixed_ind2 = find_str(admixed2);

% compute Fst between pairs of scaffold pops
% for p = 1:num_scaffold
%     for q = p+1:num_scaffold
%         par_p = false(num_nodes);
%         cur = parents(p); par_p(cur) = true;
%         while cur ~= num_nodes
%             cur = parents(cur); par_p(cur) = true;
%         end
%         cur = parents(q);
%         while ~par_p(cur)
%             cur = parents(cur);
%         end
%         fprintf('%s %s: %s, %.4f\n',node_names{p},node_names{q},node_names{cur},...
%             2*f2_all(scaffold_f2_inds(p),scaffold_f2_inds(q))/het_nodes(cur));
%     end
% end

% compare edge lengths to those inferred from heterozygosity (less accurate)
% edge_lens0 = zeros(num_edges,1);
% for edge = 1:num_edges
%     edge_lens0(edge) = (het_nodes(parents(edge))-het_nodes(edge))/2;
% end
% [edge_lens edge_lens0]

%%%%%%%%%% re-fit scaffold tree branch lengths %%%%%%%%%%

fprintf(fid,'\n');
fprintf(fid,'\t************************************\n');
fprintf(fid,'\t********** OPTIMIZE EDGES **********\n');
fprintf(fid,'\t************************************\n');

% save this to move the root back to where it should be after opt_edges
old_root1len = edge_lens(pointers(end,1));

[edge_lens,resnorm0,max_residual0] = opt_edges();
% [~,data0] = opt_edges(); use this instead of the previous line to turn off opt_edges

% correct the root branch lengths (retaining optimized sum of branches)
new_root1len = edge_lens(pointers(end,1));
new_root2len = edge_lens(pointers(end,2));
edge_lens(pointers(end,1)) = old_root1len;
edge_lens(pointers(end,2)) = new_root1len + new_root2len - old_root1len;

fit.resnorm0 = resnorm0 / COEFFS_MULT^2;
fit.max_residual0 = max_residual0 / COEFFS_MULT;

% set up variable naming for mixture equations

E1TOP_IND = num_edges + 1;
E2TOP_IND = num_edges + 2;
E3TOP_IND = num_edges + 3;
C11_IND = num_edges + 4;
C12_IND = num_edges + 5;
C2_IND = num_edges + 6;

%%%%%%%%%% optimize the first admixture %%%%%%%%%%

e1_vals = 1:(length(node_names)-1);

best_fval = 1e9;
best_e1 = 0; best_e2 = 0; best_alpha1 = 0;

%%%%% run to find best resnorm achieved by valid params %%%%%

fprintf(fid,'\n');
fprintf(fid,'\t*************************************\n');
fprintf(fid,'\t********** FIRST ADMIXTURE **********\n');
fprintf(fid,'\t*************************************\n');

fprintf(fid,'\nperforming fast search on edge');

FAST_SOLVE_WITHOUT_BOUNDS = true;

saved_fvals = zeros(length(node_names));
for e1 = e1_vals
    fprintf(fid,' %d',e1);
    for e2 = (e1+1):(length(node_names)-1)
        E1 = e1; E2 = e2;
        [alpha1,fval] = fminbnd(@fun1,0,1,optimset('display','off'));
        saved_fvals(e1,e2) = fval;
        if all(0<params & params<1) && fval < best_fval
            best_fval = fval; best_e1 = e1; best_e2 = e2; best_alpha1 = alpha1;
        end
    end
end

fprintf(fid,'\n\nbest_fval from initial fast search (without lsqlin): %g\n\n',best_fval);

FAST_SOLVE_WITHOUT_BOUNDS = false;
found_fit = false;

% edges are indexed by nodes: each node except the root has an edge to its parent
% the root is always the last branch point
for e1 = e1_vals
    fprintf(fid,'>>>>> processing edge %d (%s to its parent) <<<<<\n',e1,node_names{e1});
    for e2 = (e1+1):(length(node_names)-1)
        
        if saved_fvals(e1,e2) > best_fval * 1.1
            continue
        end
        
        E1 = e1; E2 = e2;
        [alpha1,fval] = fminbnd(@fun1,0,1,optimset('display','off'));
        if ~found_fit || fval < best_fval
            found_fit = true;
            best_fval = fval; best_e1 = e1; best_e2 = e2; best_alpha1 = alpha1;
            params1 = params; best_fval1 = best_fval;
            fun1(alpha1); % reset the coeffs
            [params,resnorm,residual] = lsqlin(coeffs(:,1:end-1),coeffs(:,end),[],[],[],[],zeros(3,1),ones(3,1),[],optimset('display','off'));
            fprintf(fid,'\n------------------------------------\n');
            fprintf(fid,'resnorm = %.5g\n\n',resnorm);
            trace_node(e1)
            trace_node(e2)
            fprintf(fid,'\n');
            fprintf(fid,'alpha1  = %.4f\n',alpha1);
            fprintf(fid,'r1      = %.4f\n',params(1));
            fprintf(fid,'r2      = %.4f\n',params(2));
            fprintf(fid,'c1      = %.4f\n',params(3));
            fprintf(fid,'\n');
            for p = 1:num_scaffold
                fprintf(fid,'%7.4f ',coeffs(p,:));
                fprintf(fid,'%16s %10.4f %10.4f\n',node_names{p},residual(p),f2_all(admixed_ind1,scaffold_f2_inds(p)));
            end
            fprintf(fid,'\n');
        end
    end
end

%%%%% swap branches if alpha < 0.5, so branch1 is the major contribution %%%%%

if best_alpha1 < 0.5
    tmp = best_e1; best_e1 = best_e2; best_e2 = tmp;
    best_alpha1 = 1-best_alpha1;
    tmp = params1(1); params1(1) = params1(2); params1(2) = tmp;
    tmp = params(1); params(1) = params(2); params(2) = tmp;
end

%%%%%%%%%% optimize the second admixture using optimal branches of first %%%%%%%%%%

if ~isempty(admixed2)

    fprintf(fid,'\n');
    fprintf(fid,'\t**************************************\n');
    fprintf(fid,'\t********** SECOND ADMIXTURE **********\n');
    fprintf(fid,'\t**************************************\n\n');

    best_fval = 1e9;
    E1 = best_e1; E2 = best_e2;

    fprintf(fid,'\nperforming fast search on edge');

    FAST_SOLVE_WITHOUT_BOUNDS = true;

    for e3 = 1:(length(node_names)-1)
        fprintf(fid,' %d',e3);
        E3 = e3;
        [alphas,fval] = fminbnd(@fun2,0,1,optimset('display','off'));
        saved_fvals(e3) = fval;
        if all(0<params & params<1) && fval < best_fval
            best_fval = fval;
        end
    end

    fprintf(fid,'\nlower bounds on resnorms achievable for e3 choices:\n');
    if fid==1
        saved_fvals(1:(length(node_names)-1))
    end

    fprintf(fid,'\n\nbest_fval from initial fast search (without lsqlin): %g\n\n',best_fval);

    FAST_SOLVE_WITHOUT_BOUNDS = false;
    found_fit = false;
    best_e3 = 0; best_alpha2 = 0;

    for e3 = 1:(length(node_names)-1)
        
        if saved_fvals(e3) > 1.1*best_fval % look at other ones that are close
            continue
        end
        
        fprintf(fid,'>>>>> processing edge %d (%s to its parent) <<<<<\n',e3,node_names{e3});
        E3 = e3;
        [alphas,fval] = fminbnd(@fun2,0,1,optimset('display','off'));
        alphas(2) = alphas(1); alphas(1) = best_alpha1;

        fprintf(fid,'                   resnorm: %.5g\n',fval);
        
        if ~found_fit || fval < best_fval
            found_fit = true;
            best_fval = fval; best_e3 = e3; best_alpha2 = alphas(2);
            display_output2(false);
        end
    end

    %%%%% fully optimize best e3 %%%%%

    if ~check_option(options,'test_2v3') % don't re-optimize when comparing resnorms of 3-way-onto-C1 and 2-way
    
    fprintf(fid,'\n');
    fprintf(fid,'\t***************************************\n');
    fprintf(fid,'\t********** OPTIMIZE ALPHA1&2 **********\n');
    fprintf(fid,'\t***************************************\n\n');

    e3 = best_e3; E3 = e3;
    init_alphas = [best_alpha1 best_alpha2];
    [alphas_orig_3way,best_fval_orig_3way] = fminsearch(@fun2,init_alphas);
    display_output2(true);
    
    %%%%% save r1, r2, r3 and re-optimize alpha1, alpha2, c11, c2 for admixed2 %%%%%
    fun2(alphas_orig_3way); % reset params to optimal
    params_orig_3way = params;
    
    if ~check_option(options,'no_opt3way')
    init_alphas = alphas_orig_3way;
    [alphas_opt,best_fval_opt] = fminsearch(@fun2_opt_3way_fracs,init_alphas);
    fun2_opt_3way_fracs(alphas_opt); % reset params to optimal
    params_opt = params;    
    end
    
    end

end

fprintf(fid,'\n');
fprintf(fid,'\t*************************************\n');
fprintf(fid,'\t********** SUMMARY RESULTS **********\n');
fprintf(fid,'\t*************************************\n\n');

fprintf(fid,'---- first admixture ----\n');
fprintf(fid,'%s\n',admixed1);
trace_node(best_e1)
trace_node(best_e2)
fprintf(fid,'%.5g\n',best_fval1/COEFFS_MULT^2);
fprintf(fid,'%.4f\n',best_alpha1);
fprintf(fid,'%.4f\n',params1(1));
fprintf(fid,'%.4f\n',params1(2));
fprintf(fid,'%.4f\n',params1(3));

fit.trace1 = get_trace(best_e1); fit.split1 = compute_split(best_e1); fit.set1 = compute_set(best_e1);
fit.trace2 = get_trace(best_e2); fit.split2 = compute_split(best_e2); fit.set2 = compute_set(best_e2);
% fit.data1 = [best_fval1/COEFFS_MULT^2; best_alpha1; params1(:)];
fit.resnorm1 = best_fval1/COEFFS_MULT^2;
fit.alpha = best_alpha1;
fit.r1 = params1(1);
fit.r2 = params1(2);
fit.c = params1(3);

fit.h1 = het_nodes(parents(best_e1));
fit.h2 = het_nodes(parents(best_e2));
fit.edge_len1 = edge_lens(best_e1);
fit.edge_len2 = edge_lens(best_e2);
fit.h_mix1 = het(find_str(admixed1));

if ~isempty(admixed2)
    
if check_option(options,'test_2v3')
fit.resnorm12 = best_fval/COEFFS_MULT^2; % for testing 2-way vs. 3-way goodness of fit, all we need is the resnorm
else

fprintf(fid,'---- second admixture ----\n');
fprintf(fid,'%s\n',admixed2);
trace_node(best_e3);
fprintf(fid,'resnorm:  %.5g\n',best_fval_orig_3way/COEFFS_MULT^2);
fprintf(fid,'alpha1:   %.4f\n',alphas_orig_3way(1));
fprintf(fid,'alpha2:   %.4f\n',alphas_orig_3way(2));
fprintf(fid,'r1:       %.4f\n',params_orig_3way(1));
fprintf(fid,'r2:       %.4f\n',params_orig_3way(2));
fprintf(fid,'r3:       %.4f\n',params_orig_3way(3));
fprintf(fid,'c11:      %.4f\n',params_orig_3way(4));
fprintf(fid,'c12:      %.4f\n',params_orig_3way(5));
fprintf(fid,'c2:       %.4f\n',params_orig_3way(6));

fit.trace3 = get_trace(best_e3); fit.split3 = compute_split(best_e3); fit.set3 = compute_set(best_e3);
% fit.data2 = [best_fval/COEFFS_MULT^2; alphas(:); params(:)];
fit.resnorm12 = best_fval_orig_3way/COEFFS_MULT^2;
fit.alpha_refit = alphas_orig_3way(1);
fit.alpha2 = alphas_orig_3way(2);
fit.r1_refit = params_orig_3way(1);
fit.r2_refit = params_orig_3way(2);
fit.r3 = params_orig_3way(3);
fit.c11 = params_orig_3way(4);
fit.c12 = params_orig_3way(5);
fit.c2 = params_orig_3way(6);

fit.h3 = het_nodes(parents(best_e3));
fit.edge_len3 = edge_lens(best_e3);
fit.h_mix2 = het(find_str(admixed2));

if ~check_option(options,'no_opt3way')
fprintf(fid,'---- optimized alpha1, alpha2, c11, c2 for admixed2 (n eqns)----\n');
fprintf(fid,'resnorm:  %.5g (%.5g pre-opt from resnorm12-resnorm1 on 3-way)\n',best_fval_opt/COEFFS_MULT^2,fit.resnorm12-fit.resnorm1);
fprintf(fid,'alpha1:   %.4f\n',alphas_opt(1));
fprintf(fid,'alpha2:   %.4f\n',alphas_opt(2));
fprintf(fid,'r1:       %.4f\n',params_opt(1));
fprintf(fid,'r2:       %.4f\n',params_opt(2));
fprintf(fid,'r3:       %.4f\n',params_opt(3));
fprintf(fid,'c11:      %.4f\n',params_opt(4));
fprintf(fid,'c12:      %.4f\n',params_opt(5));
fprintf(fid,'c2:       %.4f\n',params_opt(6));

fit.resnorm2_opt = best_fval_opt/COEFFS_MULT^2;
% replace mixture parameters with optimized versions
fit.resnorm12 = fit.resnorm1+fit.resnorm2_opt;
fit.alpha_refit = alphas_opt(1);
fit.alpha2 = alphas_opt(2);
fit.r1_refit = params_opt(1);
fit.r2_refit = params_opt(2);
fit.r3 = params_opt(3);
fit.c11 = params_opt(4);
fit.c12 = params_opt(5);
fit.c2 = params_opt(6);
end

end

end

if fid ~= 1
    fclose(fid);
end

    function option_val = check_option(options,option_name)
        default_options.print_output = 1;
        default_options.drift_units = 0;
        default_options.no_opt3way = 0;
        default_options.test_2v3 = 0;
        if isfield(options,option_name)
            option_val = options.(option_name);
        else
            option_val = default_options.(option_name);
        end
    end

    function check_admixed_pop_valid(str)
        if ~isempty(str)
            if ~any(strcmp(pop_names,str))
                error(['admixed pop ' str ' not found in pop_data.pop_names; check spelling?'])
            end
            if any(strcmp(scaffold_pop_names,str))
                error(['admixed pop ' str ' in scaffold_pop_names']);
            end
        end
    end

    function display_output2(print_eqns)
        fun2(alphas); % reset the coeffs
        if ~check_option(options,'test_2v3')
            [params,resnorm,residual] = lsqlin(coeffs(:,1:end-1),coeffs(:,end),[],[],[],[],zeros(6,1),ones(6,1),[],optimset('display','off'));
        else % leave out last equation
            Aeq = [1 0 0 0 0 0
                   0 1 0 0 0 0
                   0 0 0 1 1 0];
            beq = params1;
            warning('off','all')
            [params,resnorm,residual] = lsqlin(coeffs(:,1:end-1),coeffs(:,end),[],[],Aeq,beq,zeros(6,1),ones(6,1),[],optimset('display','off'));
            warning('on','all')
        end
        fprintf(fid,'\n------------------------------------\n');
        fprintf(fid,'resnorm = %.5g\n\n',resnorm);
        trace_node(e3)
        fprintf(fid,'\n');
        fprintf(fid,'alpha1  = %.4f\n',alphas(1));
        fprintf(fid,'alpha2  = %.4f\n',alphas(2));
        fprintf(fid,'r1      = %.4f\n',params(1));
        fprintf(fid,'r2      = %.4f\n',params(2));
        fprintf(fid,'r3      = %.4f\n',params(3));
        fprintf(fid,'c11     = %.4f\n',params(4));
        fprintf(fid,'c12     = %.4f\n',params(5));
        fprintf(fid,'c2      = %.4f\n',params(6));
        fprintf(fid,'\n');
        
        if print_eqns
            fprintf(fid,'---- equations with mix1 ----\n');
            for scaffold_ind = 1:num_scaffold
                fprintf(fid,'%7.4f ',coeffs(scaffold_ind,:));
                fprintf(fid,'%16s %10.4f %10.4f\n',node_names{scaffold_ind},residual(scaffold_ind),f2_all(admixed_ind1,scaffold_f2_inds(scaffold_ind)));
            end
            fprintf(fid,'---- equations with mix2 ----\n');
            for scaffold_ind = 1:num_scaffold
                fprintf(fid,'%7.4f ',coeffs(p+num_scaffold,:));
                fprintf(fid,'%16s %10.4f %10.4f\n',node_names{scaffold_ind},residual(scaffold_ind+num_scaffold),f2_all(admixed_ind2,scaffold_f2_inds(scaffold_ind)));
            end
            fprintf(fid,'---- equation with mix1 and mix2 ----\n');
            fprintf(fid,'%7.4f ',coeffs(end,:));
            fprintf(fid,'%16s %10.4f %10.4f\n','',residual(2*num_scaffold+1),f2_all(admixed_ind1,admixed_ind2));
            fprintf(fid,'\n');
        end
    end

    function coeffs = coeffs2root(node)
        coeffs = zeros(numel(pointers)+6,1);
        while parents(node) ~= 0
            coeffs(node) = 1;
            node = parents(node);
        end
    end

    function coeffs = coeffs2root_mix1(e1,e2,alpha1)
        coeffs = zeros(numel(pointers)+6,1);
        coeffs(C12_IND) = 1;
        coeffs(C11_IND) = 1;
        coeffs(E1TOP_IND) = alpha1;
        coeffs = coeffs + alpha1 * coeffs2root(parents(e1));
        coeffs(E2TOP_IND) = 1-alpha1;
        coeffs = coeffs + (1-alpha1) * coeffs2root(parents(e2));
    end

    function coeffs = coeffs2root_mix2(e1,e2,e3,alpha1,alpha2)
        coeffs = zeros(numel(pointers)+6,1);
        coeffs(C2_IND) = 1;
        
        % contribution of mix1, minus the c12 drift after the mix2 event
        coeffs_mix1 = coeffs2root_mix1(e1,e2,alpha1);
        coeffs_mix1(C12_IND) = 0;
        
        coeffs = coeffs + alpha2 * coeffs_mix1;
        
        coeffs(E3TOP_IND) = (1-alpha2);
        coeffs = coeffs + (1-alpha2) * coeffs2root(parents(e3));
    end

    function trace_str = get_trace(child)
        trace_str = '';
        while child > num_scaffold
            trace_str = strcat(trace_str,sprintf('%d-',child));
            child = pointers(child-num_scaffold,1); 
        end
        trace_str = strcat(trace_str,node_names{child});
    end

    function trace_node(child)
        fprintf(fid,'trace: %-20sset: %s\n',get_trace(child),compute_set(child));
    end

    function split = compute_split(edge)
        split = false(num_scaffold+1,1);
        for scaffold_pop = 1:num_scaffold
            tree_ind = find(strcmp(node_names,scaffold_pop_names{scaffold_pop}),1);
            split(scaffold_pop) = tree_dist(tree_ind,edge) < tree_dist(tree_ind,parents(edge));
        end
        split(scaffold_pop+1) = tree_dist(num_edges+1,edge) < tree_dist(num_edges+1,parents(edge)); % root
        % canonicalize the split
        if sum(~split) < sum(split) || sum(~split)==sum(split) && split(1)
            split = ~split;
        end
    end

    function set = compute_set(edge)
        split = compute_split(edge);
        set = '';
        scaffold_root_pop_names = {scaffold_pop_names{:} 'root'};
        for pop_name_cell = sort(scaffold_root_pop_names(split))
            if isempty(set)
                set = ['(' pop_name_cell{1}];
            else
                set = [set ',' pop_name_cell{1}];
            end
        end
        set = [set ')'];
        if sum(split)==1 % get rid of parentheses
            set = set(2:end-1);
        end
    end

    function resnorm = fun1(alpha1)
        coeffs1 = coeffs2root_mix1(E1,E2,alpha1);
        coeffs = zeros(num_scaffold,4);
        for scaffold_pop = 1:num_scaffold
            coeffs_scaffold = coeffs2root(scaffold_pop);
            cdiff = coeffs1 - coeffs_scaffold;
            const_term = f2_all(admixed_ind1,scaffold_f2_inds(scaffold_pop))...
                - sum(edge_lens(1:num_edges).*cdiff(1:num_edges).^2);

            % r1, r2, c11+c12, const
            coeffs(scaffold_pop,1) = edge_lens(E1) * ((cdiff(E1TOP_IND)+cdiff(E1))^2 - cdiff(E1)^2);
            coeffs(scaffold_pop,2) = edge_lens(E2) * ((cdiff(E2TOP_IND)+cdiff(E2))^2 - cdiff(E2)^2);
            coeffs(scaffold_pop,3) = cdiff(C11_IND)^2;
            coeffs(scaffold_pop,4) = const_term;
        end
        coeffs = coeffs*COEFFS_MULT;

        if FAST_SOLVE_WITHOUT_BOUNDS
            warning('off','all')
            params = coeffs(:,1:end-1)\coeffs(:,end);
            warning('on','all')
            resnorm = norm(coeffs(:,1:end-1)*params - coeffs(:,end))^2;
        else
            [params,resnorm] = lsqlin(coeffs(:,1:end-1),coeffs(:,end),[],[],[],[],zeros(3,1),ones(3,1),[],optimset('display','off'));
        end
    end

    function resnorm = fun2(alphas)
        if length(alphas) == 1
            alphas(2) = alphas(1);
            alphas(1) = best_alpha1;
        end
        coeffs1 = coeffs2root_mix1(E1,E2,alphas(1));
        coeffs2 = coeffs2root_mix2(E1,E2,E3,alphas(1),alphas(2));
        if ~check_option(options,'test_2v3')
            num_eqns = 2*num_scaffold+1;
        else
            num_eqns = 2*num_scaffold; % leave out last equation
        end
        coeffs = zeros(num_eqns,7);
        for eqno = 1:2*num_eqns
            if eqno <= num_scaffold
                scaffold_pop = eqno;
                coeffs_scaffold = coeffs2root(scaffold_pop);
                cdiff = coeffs1 - coeffs_scaffold;
                const_term = f2_all(admixed_ind1,scaffold_f2_inds(scaffold_pop));
            elseif eqno <= 2*num_scaffold
                scaffold_pop = eqno - num_scaffold;
                coeffs_scaffold = coeffs2root(scaffold_pop);
                cdiff = coeffs2 - coeffs_scaffold;
                const_term = f2_all(admixed_ind2,scaffold_f2_inds(scaffold_pop));
            else
                cdiff = coeffs1 - coeffs2;
                const_term = f2_all(admixed_ind1,admixed_ind2);
            end            
            const_term = const_term - sum(edge_lens(1:num_edges).*cdiff(1:num_edges).^2);

            % r1, r2, r3, c11, c12, c2, const
            coeffs(eqno,1) = edge_lens(E1) * ((cdiff(E1TOP_IND)+cdiff(E1))^2 - cdiff(E1)^2);
            coeffs(eqno,2) = edge_lens(E2) * ((cdiff(E2TOP_IND)+cdiff(E2))^2 - cdiff(E2)^2);
            coeffs(eqno,3) = edge_lens(E3) * ((cdiff(E3TOP_IND)+cdiff(E3))^2 - cdiff(E3)^2);
            coeffs(eqno,4) = cdiff(C11_IND)^2;
            coeffs(eqno,5) = cdiff(C12_IND)^2;
            coeffs(eqno,6) = cdiff(C2_IND)^2;
            coeffs(eqno,7) = const_term;
        end
        coeffs = coeffs*COEFFS_MULT;
        
        if FAST_SOLVE_WITHOUT_BOUNDS
            warning('off','all')
            params = coeffs(:,1:end-1)\coeffs(:,end);
            warning('on','all')
            resnorm = norm(coeffs(:,1:end-1)*params - coeffs(:,end))^2;
        else
            if ~check_option(options,'test_2v3')
                [params,resnorm] = lsqlin(coeffs(:,1:end-1),coeffs(:,end),[],[],[],[],zeros(6,1),ones(6,1),[],optimset('display','off'));
            else % leave out last equation
                Aeq = [1 0 0 0 0 0
                       0 1 0 0 0 0
                       0 0 0 1 1 0];
                beq = params1;
                warning('off','all')
                [params,resnorm] = lsqlin(coeffs(:,1:end-1),coeffs(:,end),[],[],Aeq,beq,zeros(6,1),ones(6,1),[],optimset('display','off'));
                warning('on','all')
            end
        end
    end

    % optimize alpha1, alpha2, c11, c2; hold fixed r1, r2, r3, and c12=0
    function resnorm = fun2_opt_3way_fracs(alphas)
        if length(alphas) == 1
            alphas(2) = alphas(1);
            alphas(1) = best_alpha1;
        end
        coeffs1 = coeffs2root_mix1(E1,E2,alphas(1));
        coeffs2 = coeffs2root_mix2(E1,E2,E3,alphas(1),alphas(2));
%         coeffs = zeros(2*num_scaffold+1,7);
        coeffs = zeros(num_scaffold,7); % use only the f2(C2,scaffold) eqns
        for eqno = 1:num_scaffold
                scaffold_pop = eqno;
                coeffs_scaffold = coeffs2root(scaffold_pop);
                cdiff = coeffs2 - coeffs_scaffold;
                const_term = f2_all(admixed_ind2,scaffold_f2_inds(scaffold_pop));
                
            const_term = const_term - sum(edge_lens(1:num_edges).*cdiff(1:num_edges).^2);

            % r1, r2, r3, c11, c12, c2, const
            coeffs(eqno,1) = edge_lens(E1) * ((cdiff(E1TOP_IND)+cdiff(E1))^2 - cdiff(E1)^2);
            coeffs(eqno,2) = edge_lens(E2) * ((cdiff(E2TOP_IND)+cdiff(E2))^2 - cdiff(E2)^2);
            coeffs(eqno,3) = edge_lens(E3) * ((cdiff(E3TOP_IND)+cdiff(E3))^2 - cdiff(E3)^2);
            coeffs(eqno,4) = cdiff(C11_IND)^2;
            coeffs(eqno,5) = cdiff(C12_IND)^2;
            coeffs(eqno,6) = cdiff(C2_IND)^2;
            coeffs(eqno,7) = const_term;
        end
        coeffs = coeffs*COEFFS_MULT;
        
        if FAST_SOLVE_WITHOUT_BOUNDS
            warning('off','all')
            params = coeffs(:,1:end-1)\coeffs(:,end);
            warning('on','all')
            resnorm = norm(coeffs(:,1:end-1)*params - coeffs(:,end))^2;
        else
            Aeq = [1 0 0 0 0 0
                   0 1 0 0 0 0
                   0 0 1 0 0 0
                   0 0 0 0 1 0];
            beq = [params_orig_3way(1:3); 0];
            warning('off','all')
            [params,resnorm] = lsqlin(coeffs(:,1:end-1),coeffs(:,end),[],[],Aeq,beq,zeros(6,1),ones(6,1),[],optimset('display','off'));
            warning('on','all')
        end
    end

    function [lsq_edge_lens,resnorm0,max_residual0] = opt_edges()
%        coeffs1 = coeffs2root_mix1(E1,E2,alphas(1));
%        coeffs2 = coeffs2root_mix2(E1,E2,E3,alphas(1),alphas(2));
%         coeffs = zeros(2*num_scaffold+1,7);
        num_eqns = num_scaffold*(num_scaffold-1)/2;
        coeffs = zeros(num_eqns,num_edges + 1);
            %we add 1 to allow for the constant term
        
        scaffold_pop1 = 1; scaffold_pop2 = 2;
        for eqno = 1:num_eqns
            coeffs_scaffold1 = coeffs2root(scaffold_pop1);
            coeffs_scaffold2 = coeffs2root(scaffold_pop2);
            cdiff = coeffs_scaffold1 - coeffs_scaffold2;
            const_term = f2_all(scaffold_f2_inds(scaffold_pop1),scaffold_f2_inds(scaffold_pop2));
            scaffold_pop2 = scaffold_pop2 + 1;
            if scaffold_pop2 > num_scaffold
                scaffold_pop1 = scaffold_pop1 + 1;
                scaffold_pop2 = scaffold_pop1 + 1;
            end
        
            
%             for e = 1:num_edges
%                 const_term = const_term - edge_lens(e) * cdiff(e)^2;
%             end
            for e = 1:num_edges
                coeffs(eqno,e) = cdiff(e)^2;
            end


            coeffs(eqno,end) = const_term;
        end
        coeffs = coeffs*COEFFS_MULT;
        [lsq_edge_lens,resnorm0,residuals] = lsqlin(coeffs(:,1:end-1),coeffs(:,end),[],[],[],[],zeros(num_edges,1),ones(num_edges,1),[],optimset('display','off'));
        max_residual0 = max(abs(residuals));
        
        fprintf(fid,'total resnorm: %.4f\n',resnorm0);
        fprintf(fid,'max residual: %.4f\n\n',max_residual0);
        for e = 1:num_edges
            fprintf(fid,'%s re-fitted: %.6f\n',node_names{e},lsq_edge_lens(e));
        end
    end

    function f2 = f2_from_str(sA,sB)
        f2 = f2_all(find_str(sA),find_str(sB));
    end

    function ind = find_str(s)
        ind = find(strcmp(pop_names,s),1);
    end

end

function [resnorm,het_nodes,reroot_dist] = fit_het(f2_all,het,pop_names,scaffold_pop_names,tree)

num_scaffold = length(scaffold_pop_names);
num_edges = 2*(num_scaffold-1);
num_nodes = num_edges+1;

node_names = get(tree,'nodenames');
pointers = get(tree,'pointers');
scaffold_f2_inds = zeros(num_scaffold,1);
for p = 1:num_scaffold
    scaffold_f2_inds(p) = find_str(node_names{p});
end
parents = zeros(num_nodes,1);
for ptr = 1:size(pointers,1)
    for child=1:2
        parents(pointers(ptr,child)) = ptr + num_scaffold;
    end
end

%%%%% system of equations for optimizing heterozygosities at current and ancestral nodes %%%%%
% - num_scaffold fits to h's
% - (num_scaffold choose 2) fits to (h1+h2+2*f2)/2 eqns
%
% coeffs = zeros(num_scaffold + num_scaffold*(num_scaffold-1)/2,num_nodes+1);
% eqno = 0;

% set the heterozygosities of the scaffold nodes
het_nodes = zeros(num_nodes,1);
for p = 1:num_scaffold
    het_nodes(p) = het(scaffold_f2_inds(p));
%     eqno = eqno+1;
%     coeffs(eqno,p) = 1;
%     coeffs(eqno,end) = het(scaffold_f2_inds(p));
end

is_descendant = false(num_scaffold,num_nodes);
for p = 1:num_scaffold
    cur_node = p;
    while cur_node ~= num_nodes % index of root
        is_descendant(p,cur_node) = true;
        cur_node = parents(cur_node);
    end
end

% compute best-fit heterozygosities for the ancestral nodes
ind_list = 1:num_scaffold;
resnorm = 0;
for anc = num_scaffold+1:num_nodes
%     node_names(anc)
%     node_names(ind_list(is_descendant(:,pointers(anc-num_scaffold,1))))
%     node_names(ind_list(is_descendant(:,pointers(anc-num_scaffold,2))))
    i_list = ind_list(is_descendant(:,pointers(anc-num_scaffold,1)));
    j_list = ind_list(is_descendant(:,pointers(anc-num_scaffold,2)));
    h_est = zeros(length(i_list),length(j_list));
    for i = 1:length(i_list)
        for j = 1:length(j_list)
            pi = i_list(i); pj = j_list(j);
            h_est(i,j) = (het_nodes(pi) + het_nodes(pj) + 2*f2_all(scaffold_f2_inds(pi),scaffold_f2_inds(pj)))/2;
%             eqno = eqno+1;
%             coeffs(eqno,anc) = 1;
%             coeffs(eqno,pi) = -0.5; coeffs(eqno,pj) = -0.5;
%             coeffs(eqno,end) = f2_all(scaffold_f2_inds(pi),scaffold_f2_inds(pj));
        end
    end
%     h_est
%     std(h_est(:))
    het_nodes(anc) = mean(h_est(:));
    resnorm = resnorm + sum((h_est(:)-mean(h_est(:))).^2);
end

%%%%% use heterozygosities computed by fitting current heterozygosities as well %%%%%
% coeffs(1:p,:) = coeffs(1:p,:)*.01;
% params = coeffs(:,1:end-1)\coeffs(:,end);
% [het_nodes params]
% het_nodes = params;
% resnorm = norm(coeffs(:,1:end-1)*params - coeffs(:,end))^2;

%%%%% compute edge lengths %%%%%
% edge_lens = zeros(num_edges,1);
% for edge = 1:num_edges
%     edge_lens(edge) = (het_nodes(parents(edge))-het_nodes(edge))/2;
% end

% get node closest to root: this is the original branch node
original_branch = find(select(tree,2),1);
reroot_dist = (het_nodes(parents(original_branch))-het_nodes(original_branch))/2;
if reroot_dist < 0
    reroot_dist = 0;
end

    function ind = find_str(s)
        ind = find(strcmp(pop_names,s),1);
    end

end
