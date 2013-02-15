function [trxcurr,pred] = ...
  TrackTwoFliesOneFrameOneROI(isforebb,dbkgdbb,pred,trxprev,nflies_per_roi,params)

% initialize
trxcurr = trxprev;
trxcurr.gmm_isbadprior = 0;
trxcurr.istouching = 0;
trxpriors = zeros(1,2);

%% fit connected components

cc = bwconncomp(isforebb);

cc.Area = cellfun(@numel,cc.PixelIdxList);
isbigenough = cc.Area >= params.minccarea;
ncc_bigenough = nnz(isbigenough);

%% only one fly
    
if nflies_per_roi == 1,
  
  % store in position 1
  i = 1;
  if ncc_bigenough == 1,
    % use this connected component
    j = find(cc.Area >= params.minccarea);
    idx = cc.PixelIdxList{j};
  else
    idx = cat(1,cc.PixelIdxList{:});
  end
  [y,x] = ind2sub(size(isforebb),idx);
  w = dbkgdbb(idx);
  [mu,S] = weighted_mean_cov([x(:),y(:)],w);
  [trxcurr.a(i),trxcurr.b(i),trxcurr.theta(i)] = cov2ell(S);
  trxcurr.x(i) = mu(1);
  trxcurr.y(i) = mu(2);
  trxcurr.area(i) = numel(cc.PixelIdxList{j});
  trxpriors(i) = 1;
  
  trxcurr.istouching = 0;
  order = [1,2];    

else
  
  %% two flies
  
  % fit ellipses

  % separate connected components
  if ncc_bigenough == 2,
    ccidx = find(isbigenough);
    for i = 1:2,
      j = ccidx(i);
      [y,x] = ind2sub(size(isforebb),cc.PixelIdxList{j});
      w = dbkgdbb(cc.PixelIdxList{j});
      [mu,S] = weighted_mean_cov([x(:),y(:)],w);
      [trxcurr.a(i),trxcurr.b(i),trxcurr.theta(i)] = cov2ell(S);
      trxcurr.x(i) = mu(1);
      trxcurr.y(i) = mu(2);
      trxcurr.area(i) = numel(cc.PixelIdxList{j});
      trxpriors(i) = sum(w);
    end
    trxpriors = trxpriors / sum(trxpriors);
    trxcurr.istouching = 0;

  elseif cc.NumObjects == 0,
    warning('No flies detected in ROI');
    
  else
    
    % GMM clustering of all foreground pixels
    idx = cat(1,cc.PixelIdxList{:});
    [y,x] = ind2sub(size(isforebb),idx);
    w = dbkgdbb(idx);
    if pred.isfirstframe,
      [mu,S,trxpriors,post,nll,mixprev] = mygmm([x,y],2,...
        'Replicates',params.gmmem_nrestarts_firstframe,...
        'precision',params.gmmem_precision,...
        'MaxIters',params.gmmem_maxiters,...
        'weights',w);
    else
      [mu,S,trxpriors,post,nll,mixprev] = mygmm([x,y],2,...
        'Start',pred.mix,...
        'precision',params.gmmem_precision,...
        'MaxIters',params.gmmem_maxiters,...
        'weights',w);
      
      % check that all went well
      if any(trxpriors <= params.gmmem_min_obsprior),

        fprintf('Bad prior found, trying to reinitialize\n');
        trxcurr.gmm_isbadprior = true;

        [mu1,S1,obspriors1,post1,nll1,mixprev1] = mygmm([x,y],2,...
          'Replicates',params.gmmem_nrestarts_firstframe,...
          'precision',params.gmmem_precision,...
          'MaxIters',params.gmmem_maxiters,...
          'weights',w);

        if nll1 <= nll,
          fprintf('Using results from reinitialization, which improve nll by %f\n',nll-nll1);
          mu = mu1;
          S = S1;
          trxpriors = obspriors1(:)';
          post = post1;
          nll = nll1;
          mixprev = mixprev1;
        else
          fprintf('Reinitialization does not improve nll.\n');
        end
       
      end
    end
      
    for i = 1:2,
      [trxcurr.a(i),trxcurr.b(i),trxcurr.theta(i)] = cov2ell(S(:,:,i));
      trxcurr.x(i) = mu(i,1);
      trxcurr.y(i) = mu(i,2);
      trxcurr.area(i) = sum(post(:,i));
    end
    trxcurr.istouching = true;
      
  end
   
  % match
    
  order = 1:2;
  if ~pred.isfirstframe,      

    besterr = inf;
    for i = 1:2,
      if i == 1,
        ordercurr = [1,2];
      else
        ordercurr = [2,1];
      end
      
      dpos2 = (pred.x-trxcurr.x(ordercurr)).^2 + (pred.y-trxcurr.y(ordercurr)).^2;
      dtheta = abs(modrange(pred.theta-trxcurr.theta(ordercurr),-pi/2,pi/2));
      darea = abs(pred.area-trxcurr.area(ordercurr));
        
      errcurr = sqrt(sum(dpos2))*params.err_weightpos + ...
        sqrt(sum(dtheta.^2))*params.err_weighttheta + ...
        sqrt(sum(darea.^2))*params.err_weightarea;
        
      if errcurr < besterr,
        order = ordercurr;
        besterr = errcurr;
      end
    end
      
  end
 
  trxcurr.x = trxcurr.x(order);
  trxcurr.y = trxcurr.y(order);
  trxcurr.a = trxcurr.a(order);
  trxcurr.b = trxcurr.b(order);
  trxcurr.theta = trxcurr.theta(order);
  trxcurr.area = trxcurr.area(order);
  
end
        
%% predicted position for next frame
pred.area = trxcurr.area;
if pred.isfirstframe,
  pred.x = trxcurr.x;
  pred.y = trxcurr.y;
  pred.theta = trxcurr.theta;
else
  pred.x = (2-params.err_dampen_pos)*trxcurr.x - (1-params.err_dampen_pos)*trxprev.x;
  pred.y = (2-params.err_dampen_pos)*trxcurr.y - (1-params.err_dampen_pos)*trxprev.y;
  dtheta = modrange(trxcurr.theta-trxprev.theta,-pi/2,pi/2);
  pred.theta = trxcurr.theta+(1-params.err_dampen_theta)*dtheta;
end
   
% switch around priors, move toward .5, .5
pred.mix.priors = (1-params.err_dampen_priors)*trxpriors(order) + params.err_dampen_priors*.5;
% set centres, covars to predicted positions
pred.mix.centres = [pred.x,pred.y];
pred.mix.covars = axes2cov(trxcurr.a,trxcurr.b,pred.theta);

pred.isfirstframe = false;
