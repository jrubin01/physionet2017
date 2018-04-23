#!/usr/bin/python -tt

# # PhysioNet Test Entry
# # This script will be invoked by next.sh
# ### Input: an ECG record to process (.mat). 
# ### Output: appends classification result to answers.txt

import os
import sys
import scipy.io
import commands
import time
import datetime
import pickle
import numpy as np
import pandas as pd
from scipy import signal
from StringIO import StringIO
from datetime import timedelta

np.random.seed()

from keras.models import Sequential, model_from_json
from keras.optimizers import adam

# Load and start Matlab Engine

import matlab.engine
eng = matlab.engine.start_matlab()
eng.addpath(r'matlab/', nargout=0)
eng.addpath(r'matlab/no_classifier_with_qrs/', nargout=0)

no_threshold1 = 0.3
no_threshold2 = 0.7

do_normalize = True
do_reshape = False

model_name = os.path.join('cnn_15', '2017-08-13_09.45.48.678287')
weights_version = 'weights-08-0.88.hdf5'

backup_model_name = os.path.join('cnn_9', '2017-07-15_19.24.02.203737')
backup_weights_version = 'weights-06-0.90.hdf5'

segment_size = 15.0
backup_segment_size = 9.0
Fs = 300.0
spec_threshold = 30

# Load saved model

def load_model(model_name, weights_version):
    model_path = os.path.join('saved_models', model_name)
    model = model_from_json(open(os.path.join(model_path, 'model.json')).read())
    model.load_weights(os.path.join(model_path, weights_version))
    model.compile(loss='categorical_crossentropy', optimizer=adam(0.001), metrics=['accuracy'])
    return model
        
def load_scaler(model_name):
    model_path = os.path.join('saved_models', model_name)
    scaler_deep_mean = np.load(os.path.join(model_path, 'scaler_deep_mean.npy'))
    scaler_deep_std = np.load(os.path.join(model_path, 'scaler_deep_std.npy'))
    return scaler_deep_mean, scaler_deep_std
    
def get_prediction(index):
    if index == 0:
        return 'A'
    if index == 1:
        return 'N'
    if index == 2:
        return 'O'
    if index == 3:
        return '~'

def add_seconds(start, secs):
    '''
    Helper function to add the specified number of seconds to a start time.
    Requires conversion to full datetime format before adding timedelta and conversion back to time
    '''
    full_datetime = datetime.datetime(100, 1, 1, start.hour, start.minute, start.second, start.microsecond)
    return (full_datetime + timedelta(seconds=secs)).time()
         
def make_prediction(datapath):
    sqi = eng.sqi_calculator_modified('{}.mat'.format(datapath), nargout=1)
    
    # Load the model
    model = load_model(model_name, weights_version)
    backup_model = load_model(backup_model_name, backup_weights_version)
    
    # Load scaler
    if do_normalize:
        scaler_deep_mean, scaler_deep_std = load_scaler(model_name)
        backup_scaler_deep_mean, backup_scaler_deep_std = load_scaler(backup_model_name)
        
    if sqi == 0.0:
        pred = '~'
        print pred, 'sqi'
    else:
        # Load ECG data
        mat = scipy.io.loadmat(datapath)
        ecg = np.reshape(mat['val'], len(mat['val'][0]))
        length = len(ecg) / Fs
        
        # Perform QRS detection
        cmd = 'gqrs -r {}'.format(datapath)
        (status, output) = commands.getstatusoutput(cmd)
        cmd = 'rdann -r {} -a qrs -e -v'.format(datapath)
        (status, output) = commands.getstatusoutput(cmd)
        qrs = pd.read_csv(StringIO(output), delim_whitespace=True)
        
        if not qrs.empty:
            # Construct spectrogram
            f, t, Sxx = signal.spectrogram(ecg, fs=Fs, nperseg=75, nfft=128, noverlap=64)
            time_window = np.round(t[1] - t[0], 2)
            delta = int(segment_size / time_window)
            backup_delta = int(backup_segment_size / time_window)
            
            index = 0
            X_deep = []
            X_deep_backup = []
            for start in qrs.Elapsed.values:
                start = datetime.datetime.strptime(start, '%M:%S.%f').time()
                end = add_seconds(start, segment_size)
                end_time = (end.minute * 60) + end.second + (end.microsecond / 1e6)
                backup_end = add_seconds(start, backup_segment_size)
                backup_end_time = (backup_end.minute * 60) + backup_end.second + (backup_end.microsecond / 1e6)
        
                # Check segment end time is within the length of the recording
                if end_time < length:
                    ########################
                    # Spectrogram segemtns #
                    ########################
            
                    for indx in range(len(t)-1):
                        # Hack: ignoring anything over 60 seconds, as it will cause error
                        # with time conversion. No record exceeds 61 seconds, so this is okay for now
                        if t[indx+1] < 60.0:
                            lower = datetime.datetime.strptime(str(round(t[indx], 3)), '%S.%f').time()
                            upper = datetime.datetime.strptime(str(round(t[indx+1], 3)), '%S.%f').time()

                            if start > lower and start < upper:
                                # Only consider complete segements
                                if Sxx[:,indx:indx+delta].shape[1] == delta:
                                    X_deep.append(np.log(Sxx[:20,indx:indx+delta]))
                                if Sxx[:,indx:indx+backup_delta].shape[1] == backup_delta:
                                    X_deep_backup.append(np.log(Sxx[:20,indx:indx+backup_delta]))
                                break
                elif backup_end_time < length:
                    # Remaining backup segments
            
                    ########################
                    # Spectrogram segemtns #
                    ########################
            
                    for indx in range(len(t)-1):
                        # Hack: ignoring anything over 60 seconds, as it will cause error
                        # with time conversion. No record exceeds 61 seconds, so this is okay for now
                        if t[indx+1] < 60.0:
                            lower = datetime.datetime.strptime(str(round(t[indx], 3)), '%S.%f').time()
                            upper = datetime.datetime.strptime(str(round(t[indx+1], 3)), '%S.%f').time()

                            if start > lower and start < upper:
                                # Only consider complete segements
                                if Sxx[:,indx:indx+backup_delta].shape[1] == backup_delta:
                                    X_deep_backup.append(np.log(Sxx[:20,indx:indx+backup_delta]))
                                    break
            
            if len(X_deep) > 0:
                # Limit spectrograms if too many
                if len(X_deep) > spec_threshold:
                    X_deep = X_deep[::3]
                
                if do_normalize:
                    for i in range(len(X_deep)):
                        X_deep[i] = (X_deep[i] - scaler_deep_mean) / scaler_deep_std
                
                # Reshape to 4-channel input
                if do_reshape:
                    X_multi = []
                    for i in range(len(X_deep)):
                        X_multi.append(np.array([X_deep[i][0:5], X_deep[i][5:10], X_deep[i][10:15], X_deep[i][15:20]]))
                    X_deep = np.array(X_multi)
                else:
                    X_deep = np.expand_dims(X_deep, 1)
            
                predictions = model.predict(X_deep)
                probs = np.mean(predictions, axis=0)
                pred = get_prediction(np.argmax(probs))
                print pred, 'main'
            else:
                # Not enough segments, use the backup model
                if len(X_deep_backup) > 0:
                
                    # Limit spectrograms if too many
                    if len(X_deep_backup) > spec_threshold:
                        X_deep_backup = X_deep_backup[::3]
                                
                    if do_normalize:
                        for i in range(len(X_deep_backup)):
                            X_deep_backup[i] = (X_deep_backup[i] - backup_scaler_deep_mean) / backup_scaler_deep_std
                    X_deep_backup = np.expand_dims(X_deep_backup, 1)
            
                    predictions = backup_model.predict(X_deep_backup)
                    probs = np.mean(predictions, axis=0)
                    pred = get_prediction(np.argmax(probs))
                    print pred, 'backup'
                else:
                    pred = '~'
                    print pred, 'none'
        else:
            pred = '~'
            print pred, 'none'
    	
	# NO classifier
    if pred == 'N' or pred == 'O':
        if abs(probs[2] - probs[1]) <= no_threshold1:
            no_pred = eng.classifierNO(datapath, nargout=1)            
            if no_pred >= no_threshold2:
                pred = 'N'
                print pred, 'no_classifier'
            else:
                pred = 'O'
                print pred, 'no_classifier'
                
	# Write to answer.txt
    f = open('answers.txt', 'a')
    filename = os.path.basename(datapath)
    f.write('{},{}\n'.format(filename, pred))
    f.close()
	
                
def main():
	if len(sys.argv) != 2:
		print 'usage: ./predict.py datapath'
		sys.exit(1)
	
	datapath = sys.argv[1]
	make_prediction(datapath)
	eng.quit()
	
if __name__ == '__main__':
	main()

