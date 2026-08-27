clc; clear; close all;

%% ===== User settings =====
infile = 'harvard_clean.wav'; % set to your file (or leave '' to synthesize)
targetSNR_dB = 5; % SNR of noisy signal (try 0, 5, 10)

%% ===== Read input (or synthesize if file missing) =====
if ~isempty(infile) && exist(infile,'file')
    [x, fs] = audioread(infile);
    if size(x,2) > 1
        x = mean(x,2);
    end
else
    fs = 16000;
    t = (0:1/fs:5).';
    x = 0.9*sin(2*pi*120*t) .* ...
        (0.5 + 0.5*sin(2*pi*0.25*t));
    x = x + 0.4*sin(2*pi*(300 + ...
        50*sin(2*pi*0.1*t)).*t);
end

x = x(:);
x = x / (max(abs(x))+eps);
N = numel(x);

%% ===== Add white Gaussian noise to reach target SNR =====
rng(0); % reproducible
n = randn(size(x));

Px = mean(x.^2);
Pn = mean(n.^2);

scale = sqrt(Px/(Pn*10^(targetSNR_dB/10)));
n = scale * n;

noisy = x + n;
noisy = noisy / (max(abs(noisy))+eps);

%% ===== STFT parameters =====
winLen = 1024;
hop = winLen/4; % 75% overlap
win = hamming(winLen,'periodic');

%% ===== Compute STFTs (toolbox-free) =====
Sx = stft_custom(x, winLen, hop, win);
Sn = stft_custom(n, winLen, hop, win);
Sy = stft_custom(noisy, winLen, hop, win);

%% ===== Noise PSD estimate =====
noiseSpecPow = mean(abs(Sn).^2, 2);

%% ===== 1) Frequency-domain Wiener filter =====
Yw = zeros(size(Sy));
Gfloor = 0.01;

for k = 1:size(Sy,2)

    Y2 = abs(Sy(:,k)).^2;

    % A-posteriori SNR
    gamma = Y2 ./ (noiseSpecPow + eps);

    if k == 1

        % Initial a-priori SNR
        xi = max(gamma - 1, 0);

    else

        % Decision-directed update
        prev_spec = abs(Yw(:,k-1)).^2;

        xi = 0.98 * ...
            (prev_spec./(noiseSpecPow+eps)) ...
            + 0.02 * max(gamma-1,0);

    end

    % Wiener gain
    H = xi ./ (1 + xi);

    % Gain floor
    H = max(H, Gfloor);

    Yw(:,k) = H .* Sy(:,k);

end

y_wiener = istft_custom(Yw, winLen, hop, win, N);
y_wiener = y_wiener / (max(abs(y_wiener))+eps);

%% ===== 2) Spectral Gating =====
alpha = 4.0; % over-subtraction
beta = 0.01; % gain floor

Ysg = zeros(size(Sy));

for k = 1:size(Sy,2)

    Y2 = abs(Sy(:,k)).^2;

    % Soft gate
    mask = max( ...
        (Y2 - alpha*noiseSpecPow) ./ ...
        (Y2 + eps), beta);

    mask = min(mask, 1);

    Ysg(:,k) = mask .* Sy(:,k);

end

y_specgate = istft_custom(Ysg, winLen, hop, win, N);
y_specgate = y_specgate / (max(abs(y_specgate))+eps);

%% ===== Save audio =====
audiowrite('noisy.wav', noisy, fs);
audiowrite('denoised_wiener.wav', y_wiener, fs);
audiowrite('denoised_sg.wav', y_specgate, fs);

%% ===== Plot waveforms =====
t = (0:N-1)/fs;

figure('Name','Waveforms', ...
       'Color','w', ...
       'Position',[100 100 800 700]);

subplot(4,1,1);
plot(t,x);
title('Clean (original)');
xlim([0 t(end)]);
ylabel('Amp');

subplot(4,1,2);
plot(t,noisy);
title(sprintf('Noisy (SNR = %g dB)',targetSNR_dB));
xlim([0 t(end)]);
ylabel('Amp');

subplot(4,1,3);
plot(t,y_wiener);
title('Denoised - Wiener');
xlim([0 t(end)]);
ylabel('Amp');

subplot(4,1,4);
plot(t,y_specgate);
title('Denoised - Spectral Gating');
xlim([0 t(end)]);
xlabel('Time (s)');
ylabel('Amp');

%% ===== Spectrograms =====
figure('Name','Spectrograms', ...
       'Color','w', ...
       'Position',[150 150 1000 600]);

subplot(2,2,1);
spectrogram(x,winLen,winLen-hop,winLen,fs,'yaxis');
title('Clean');

subplot(2,2,2);
spectrogram(noisy,winLen,winLen-hop,winLen,fs,'yaxis');
title('Noisy');

subplot(2,2,3);
spectrogram(y_wiener,winLen,winLen-hop,winLen,fs,'yaxis');
title('Wiener');

subplot(2,2,4);
spectrogram(y_specgate,winLen,winLen-hop,winLen,fs,'yaxis');
title('Spectral Gating');

%% ===== Done =====
fprintf(['WAVs written: noisy.wav, ' ...
         'denoised_wiener.wav, denoised_sg.wav\n']);

%% ===== Helper Function: STFT =====
function S = stft_custom(x, frame_len, hop, win)

    % Returns one-sided STFT
    x = x(:);

    L = length(x);

    nframes = ceil((L-frame_len)/hop) + 1;

    pad = (nframes-1)*hop + frame_len - L;

    xpad = [x; zeros(pad,1)];

    nfft = frame_len;

    nk = nfft/2 + 1;

    S = zeros(nk,nframes);

    for i = 1:nframes

        idx = (1:frame_len) + (i-1)*hop;

        frame = xpad(idx) .* win;

        Xf = fft(frame,nfft);

        S(:,i) = Xf(1:nk);

    end

end

%% ===== Helper Function: ISTFT =====
function x = istft_custom(S, frame_len, hop, win, target_len)

    % Inverse of stft_custom

    nk = size(S,1);
    nframes = size(S,2);

    nfft = frame_len;

    % Reconstruct full spectrum
    Sfull = zeros(nfft,nframes);

    Sfull(1:nk,:) = S;

    if mod(nfft,2) == 0

        Sfull(nk+1:end,:) = ...
            conj(flipud(S(2:nk-1,:)));

    else

        Sfull(nk+1:end,:) = ...
            conj(flipud(S(2:nk,:)));

    end

    xlen = (nframes-1)*hop + frame_len;

    xrec = zeros(xlen,1);

    win_sum = zeros(xlen,1);

    for i = 1:nframes

        frame_spec = Sfull(:,i);

        frame = real(ifft(frame_spec,nfft));

        frame = frame(1:frame_len) .* win;

        idx = (1:frame_len) + (i-1)*hop;

        xrec(idx) = xrec(idx) + frame;

        win_sum(idx) = win_sum(idx) + win.^2;

    end

    nz = win_sum > 1e-8;

    xrec(nz) = xrec(nz) ./ win_sum(nz);

    x = xrec(1:min(target_len,length(xrec)));

end
