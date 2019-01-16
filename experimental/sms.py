# we import the Twilio client from the dependency we just installed
from twilio.rest import Client

# the following line needs your Twilio Account SID and Auth Token
client = Client("ACf0fc445a33885560fd884cc8a245c7a1", "bba19a693f1a1bca2bd9f6ddb3031d31")

# change the "from_" number to your Twilio number and the "to" number
# to the phone number you signed up for Twilio with, or upgrade your
# account to send SMS to any phone number
client.messages.create(to="+16048180169", 
                       from_="+16042391598", 
                       body="Hello from Python twilio testing!")
