import os
import sys
import time
import numpy as np
import pickle as pkl
import tensorflow.compat.v1 as tf
#import tensorflow as tf
from utils import *
from sankey import *
from calculatePercent import *
from tensorflow.python.saved_model import tag_constants
from models import scGCN
#sys.stdout = open("output_log.txt", "w")
import warnings
warnings.filterwarnings("ignore")
# Set random seed
seed = 123
np.random.seed(seed)
tf.compat.v1.set_random_seed(seed)
tf.set_random_seed(seed)
# Settings
flags = tf.app.flags
FLAGS = flags.FLAGS
flags.DEFINE_string('dataset', 'input', 'data dir')
flags.DEFINE_string('output', 'results', 'predicted results')
flags.DEFINE_bool('graph', True, 'select the optional graph.')
flags.DEFINE_string('model', 'scGCN','Model string.') 
flags.DEFINE_float('learning_rate', 0.01, 'Initial learning rate.')
flags.DEFINE_integer('epochs', 200, 'Number of epochs to train.')
flags.DEFINE_integer('hidden1', 32, 'Number of units in hidden layer 1.')
#flags.DEFINE_integer('hidden2', 32, 'Number of units in hidden layer 2.')
flags.DEFINE_float('dropout', 0, 'Dropout rate (1 - keep probability).')
flags.DEFINE_float('weight_decay', 0,
                   'Weight for L2 loss on embedding matrix.')
flags.DEFINE_integer('early_stopping', 10,
                     'Tolerance for early stopping (# of epochs).')
flags.DEFINE_integer('max_degree', 3, 'Maximum Chebyshev polynomial degree.')
# Load data
final_results=pd.DataFrame()
final_results["original"]=[]
final_results["predicted"]=[]
original_list=[]
predicted_list=[]
files= os.listdir(FLAGS.dataset)
'''
##测试多线程
import threading
from multiprocessing.pool import ThreadPool
threads=[]
class myThread(threading.Thread):
        def __init__(self,func,args):
                Thread.__init__(self)
                self.func=func
                self.args=args
        def run(self):
                self.func(*self.args)
for file in files:
        threads.append(myThread(processSCGCN,(file)))
        threads[-1].start()
for thread in threads:
        thread.join()
'''
def processSCGCN(qur_file):
#for file in files:
	file=str(qur_file)
	print(file+"processing>>>")
	position=FLAGS.dataset+'/'+file
	adj, features, labels_binary_train, labels_binary_val, labels_binary_test, train_mask, pred_mask, val_mask, test_mask, new_label, true_label, index_guide,rename = load_data(position,rgraph=FLAGS.graph)
	support = [preprocess_adj(adj)]
	num_supports = 1
	model_func = scGCN

# Define placeholders
	tf.compat.v1.disable_eager_execution()
	placeholders = {
		'support':
		[tf.sparse_placeholder(tf.float32) for _ in range(num_supports)],
		'features':
		tf.sparse_placeholder(tf.float32,
			shape=tf.constant(features[2], dtype=tf.int64)),
		'labels':
			tf.placeholder(tf.float32, shape=(None, labels_binary_train.shape[1])),
		'labels_mask':
			tf.placeholder(tf.int32),
		'dropout':
			tf.placeholder_with_default(0., shape=()),
		'num_features_nonzero':
			tf.placeholder(tf.int32)  # helper variable for sparse dropout
	}
# Create model
	model = model_func(placeholders, input_dim=features[2][1], logging=True)


# Define model evaluation function
	def evaluate(features, support, labels, mask, placeholders):
		t_test = time.time()
		feed_dict_val = construct_feed_dict(features, support, labels, mask,
			placeholders)
		outs_val = sess.run([model.loss, model.accuracy], feed_dict=feed_dict_val)
		return outs_val[0], outs_val[1], (time.time() - t_test)


# Initialize session
	sess = tf.Session()
# Init variables
	sess.run(tf.global_variables_initializer())
	train_accuracy = []
	train_loss = []
	val_accuracy = []
	val_loss = []
	test_accuracy = []
	test_loss = []

# Train model

#configurate checkpoint directory to save intermediate model training weights
	saver = tf.train.Saver()
	save_dir = 'checkpoints/'+file+"/"
	if not os.path.exists(save_dir):
		os.makedirs(save_dir)

	save_path = os.path.join(save_dir, 'best_validation')

	for epoch in range(FLAGS.epochs):
		t = time.time()
    # Construct feed dictionary
		feed_dict = construct_feed_dict(features, support, labels_binary_train,
                                    train_mask, placeholders)
		feed_dict.update({placeholders['dropout']: FLAGS.dropout})
    # Training step
		outs = sess.run([model.opt_op, model.loss, model.accuracy],
                    feed_dict=feed_dict)
		train_accuracy.append(outs[2])
		train_loss.append(outs[1])
    # Validation
		cost, acc, duration = evaluate(features, support, labels_binary_val,
                                   val_mask, placeholders)
		val_loss.append(cost)
		val_accuracy.append(acc)
		test_cost, test_acc, test_duration = evaluate(features, support,
                                                  labels_binary_test,
                                                  test_mask, placeholders)
		test_accuracy.append(test_acc)
		test_loss.append(test_cost)
		saver.save(sess=sess, save_path=save_path)
		print("Epoch:", '%04d' % (epoch + 1), "train_loss=",
			"{:.5f}".format(outs[1]), "train_acc=", "{:.5f}".format(outs[2]),
			"val_loss=", "{:.5f}".format(cost), "val_acc=", "{:.5f}".format(acc),
			"time=", "{:.5f}".format(time.time() - t))
		if epoch > FLAGS.early_stopping and val_loss[-1] > np.mean(
			val_loss[-(FLAGS.early_stopping + 1):-1]):
			print("Early stopping...")
			break

	print("Finished Training....")

#'  outputs 
	all_mask = np.array([True] * len(train_mask))
	labels_binary_all = new_label

	feed_dict_all = construct_feed_dict(features, support, labels_binary_all,
                                    all_mask, placeholders)
	feed_dict_all.update({placeholders['dropout']: FLAGS.dropout})

	activation_output = sess.run(model.activations, feed_dict=feed_dict_all)[1]
	predict_output = sess.run(model.outputs, feed_dict=feed_dict_all)

#' accuracy on all masks
	ab = sess.run(tf.nn.softmax(predict_output))
	all_prediction = sess.run(
		tf.equal(sess.run(tf.argmax(ab, 1)),
			sess.run(tf.argmax(labels_binary_all, 1))))

#' accuracy on prediction masks 
	acc_train = np.sum(all_prediction[train_mask]) / np.sum(train_mask)
	acc_test = np.sum(all_prediction[test_mask]) / np.sum(test_mask)
	acc_val = np.sum(all_prediction[val_mask]) / np.sum(val_mask)
	acc_pred = np.sum(all_prediction[pred_mask]) / np.sum(pred_mask)
	print('Checking train/test/val set accuracy: {}, {}, {}'.format(
		acc_train, acc_test, acc_val))
	print('Checking pred set accuracy: {}'.format(acc_pred))
	from sklearn.metrics import accuracy_score, f1_score, precision_score, recall_score
	y_pred_train = sess.run(tf.argmax(ab, 1))[train_mask]
	y_real_train = sess.run(tf.argmax(labels_binary_all, 1))[train_mask]
	y_pred_test = sess.run(tf.argmax(ab, 1))[test_mask]
	y_real_test = sess.run(tf.argmax(labels_binary_all, 1))[test_mask]
	y_pred_val = sess.run(tf.argmax(ab, 1))[val_mask]
	y_real_val = sess.run(tf.argmax(labels_binary_all, 1))[val_mask]
        #Prediction_accuracy=np.mean([accuracy_score(y_real_train,y_pred_train,average='weighted'),accuracy_score(y_real_test,y_pred_test,average='weighted'),accuracy_score(y_real_val,y_pred_val,average='weighted')])
	F1_score=np.mean([f1_score(y_real_train, y_pred_train, average='weighted'),f1_score(y_real_test, y_pred_test, average='weighted'),f1_score(y_real_val, y_pred_val, average='weighted')])
	Precision_score=np.mean([precision_score(y_real_train, y_pred_train, average='weighted'),precision_score(y_real_test, y_pred_test, average='weighted'),precision_score(y_real_val, y_pred_val, average='weighted')])
	Recall_score=np.mean([recall_score(y_real_train, y_pred_train, average='weighted'),recall_score(y_real_test, y_pred_test, average='weighted'),recall_score(y_real_val, y_pred_val, average='weighted')])
        #print("Prediction accuracy :{}".format(acc_pred))
        #print("Prediction accuracy :{}".format(Prediction_accuracy))
	print("F1 score :{}".format(F1_score))
	print("Precision score :{}".format(Precision_score))
	print("Recall score :{}".format(Recall_score))

#' save the predicted labels of query data
	if not os.path.exists(FLAGS.output):
		os.makedirs(FLAGS.output)
	out_dir=FLAGS.output+"/"+file
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)
	scGCN_all_labels = true_label.values.flatten()  #' ground truth
	np.savetxt(out_dir + '/scGCN_all_input_labels.csv',scGCN_all_labels,delimiter=',',comments='',fmt='%s')
	np.savetxt(out_dir+'/scGCN_query_mask.csv',pred_mask,delimiter=',',comments='',fmt='%s')           
	ab = sess.run(tf.nn.softmax(predict_output))
	all_binary_prediction = sess.run(tf.argmax(ab, 1))  #' predict catogrized labels
	all_binary_labels = sess.run(tf.argmax(labels_binary_all, 1))  #' true catogrized labels
	np.savetxt(out_dir+'/scGCN_all_binary_predicted_labels.csv',all_binary_prediction,delimiter=',',comments='',fmt='%f')
	np.savetxt(out_dir+'/scGCN_index_guide.csv',index_guide,delimiter=',',comments='',fmt='%f')
	np.savetxt(out_dir+'/scGCN_all_binary_input_labels.csv',all_binary_labels,delimiter=',',comments='',fmt='%f')
# new process
	scGCN_data2_labels = scGCN_all_labels[np.where(index_guide<0)]
	data2_binary_prediction = all_binary_prediction[np.where(index_guide<0)]
	new_rename = dict(zip(rename.values(),rename.keys()))
#更换字典rename中的key与value
	data2_binary_prediction_tab=pd.DataFrame(data2_binary_prediction)
	data2_labels_prediction = data2_binary_prediction_tab.replace(new_rename)
	np.savetxt(out_dir+'/scGCN_data2_labels_2_prediction.csv',np.column_stack((scGCN_data2_labels,data2_labels_prediction )),delimiter=',',comments='',fmt='%s %s')
	threadLock.acquire()
	original_list.extend(scGCN_data2_labels)
	predicted_list.extend(data2_labels_prediction[0])
	threadLock.release()
import threading
from threading import Thread
from multiprocessing.pool import ThreadPool
threadLock=threading.Lock()
threads=[]
class myThread(threading.Thread):
        def __init__(self,func,args):
                Thread.__init__(self)
                self.func=func
                self.args=args
        def run(self):
                self.func(self.args)
for file in files:
        threads.append(myThread(processSCGCN,(file)))
        threads[-1].start()
for thread in threads:
        thread.join()
final_results["original"]=original_list
final_results["predicted"]=predicted_list
createSankey(final_results)
calculatePercent(final_results)
np.savetxt(FLAGS.output+'/final_original2predicte.csv',np.column_stack((original_list,predicted_list)),delimiter=',',comments='',fmt='%s %s')
