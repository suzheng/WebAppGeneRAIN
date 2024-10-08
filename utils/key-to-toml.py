import toml

output_file = "/Users/suzheng/Downloads/secrets.toml"
# following https://blog.streamlit.io/streamlit-firestore-continued/#part-4-securely-deploying-on-streamlit-sharing
# downloaded from https://console.firebase.google.com/project/generain-vec/settings/serviceaccounts/adminsdk
with open("/Users/suzheng/Downloads/generain-vec-firebase-adminsdk-qj5xf-cb4f820cd1.json") as json_file:
    json_text = json_file.read()

config = {"textkey": json_text}
toml_config = toml.dumps(config)

with open(output_file, "w") as target:
    target.write(toml_config)