import glob
import random
import numpy as np

import wave
import librosa
import librosa.display


def get_and_shuffle_filenames(dir_name):
    filenames = glob.glob(str(dir_name) + "/*")
    random.shuffle(filenames)
    return filenames

def get_max_length(X_unfiltered):
    X_lengths = [audio.shape[0] for _, audio in X_unfiltered]

    max_length_1 = int(np.mean(X_lengths) + 2 * np.std(X_lengths))
    max_length = int(np.floor(max_length_1 / 256) * 256)
    return max_length

def get_max_length2(X_unfiltered, std):
    X_lengths = [audio.shape[0] for _, audio in X_unfiltered]

    max_length_1 = int(np.mean(X_lengths) + std * np.std(X_lengths))
    max_length = int(np.floor(max_length_1 / 256) * 256)
    return max_length

# https://www.tensorflow.org/tutorials/audio/simple_audio
def decode_audio(file_path):
    # read file to get buffer                                                                                               
    ifile = wave.open(file_path)
    samples = ifile.getnframes()
    audio = ifile.readframes(samples)

    # convert buffer to float32 using NumPy                                                                                 
    audio_as_np_int16 = np.frombuffer(audio, dtype=np.int16)
    audio_as_np_float32 = audio_as_np_int16.astype(np.float32)
    
    # get largest absolute value
    max_val = np.max(
        np.absolute(
            [np.max(audio_as_np_float32), np.min(audio_as_np_float32)]))
    audio_normalized = audio_as_np_float32 / max_val

    return audio_normalized

def get_label(file_path):
    # label is in the filename
    parts = file_path.split("/")
    label = int(parts[2].split("_")[0])

    return label

def spect(signal, max_length):
    range_val = int(max_length/256)

    spectogram = []
    for i in range(range_val):
        window_fft = np.fft.rfft(signal[i * 256: (i + 1) * 256])[:-1]
        window_fft = list(np.abs(window_fft))
        spectogram.append(window_fft)
    spectogram = np.array(spectogram)
    spectogram = librosa.amplitude_to_db(spectogram, ref=np.max)
    spectogram = spectogram / 20
    return spectogram

def normalize_arr(array):
    mean = np.mean(array)
    std = np.std(array)
    array = array - mean
    array = array / std
    return array

def split_full(X_full, y_full):
    num_samples = len(X_full)

    tenth = int(num_samples * 0.1)
    eightyth = tenth * 8

    X_train = X_full[:eightyth]
    y_train = y_full[:eightyth]

    X_val = X_full[eightyth: eightyth + tenth]
    y_val = y_full[eightyth: eightyth + tenth]

    X_test = X_full[eightyth + tenth:]
    y_test = y_full[eightyth + tenth:]

    return [(X_train, y_train), (X_val, y_val), (X_test, y_test)]
