import json
from os.path import join, dirname
from watson_developer_cloud import SpeechToTextV1

speech_to_text = SpeechToTextV1(
	username="5e08ce2c-705b-4dec-a2ec-90549582a835",
	password="Gkek4uUEhPn7",
	x_watson_learning_opt_out=False
)

print(json.dumps(speech_to_text.models(), indent=2))