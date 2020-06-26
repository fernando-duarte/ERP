function [publag] = computePublag_(obsdate, realtime, value)

  % Diff along 2nd dimension of ~isnan; will have a 1 in the i,j entry
  % if observation i went from NaN to some value in vintage j.
  %
  % Note this is NOT exactly a release. That's because sometimes, a
  % series will be filled in retroactively, way back in its history.
  % Then, there will be lots of NaNs that turn to values, which aren't
  % new releases. That's why we need the next indicator below
  nan_to_value = diff(~isnan(value), 1, 2);
  nan_to_value = (nan_to_value == 1); % Remove -1 from when series goes from value to NaN

  % A 1 in the i,j location indicates that observation i-1 in vintage
  % j-1 was not nan.
  %
  % This helps us distinguish genuine new releases from the
  % backfilling detailed above
  lastabove_was_notnan = [zeros(1, size(value,2)-1); ...
                          ~isnan(value(1:end-1, 1:end-1))];

  rls_inds = (nan_to_value .* lastabove_was_notnan);
  rls_inds = [zeros(size(rls_inds,1),1), rls_inds];
  [rls_obs, rls_vint] = find(rls_inds);

  publag.date     = obsdate(rls_obs);
  publag.realtime = realtime(rls_vint);
  publag.mean     = mean(publag.realtime - publag.date);
  publag.median   = median(publag.realtime - publag.date);

end
