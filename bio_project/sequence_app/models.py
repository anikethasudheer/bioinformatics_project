# sequence_app/models.py
from django.db import models
from django.contrib.auth.models import User

class DNASequence(models.Model):
    user = models.ForeignKey(User, on_delete=models.CASCADE)
    sequence = models.TextField(max_length=5000) 
    length = models.PositiveIntegerField(default=0)
    gc_content = models.FloatField(default=0.0)
    submission_date = models.DateTimeField(auto_now_add=True)

    def __str__(self):
       
        return f"ID {self.id} - Length {self.length} bp"

from django.db import models
from django.contrib.auth.models import User

class SavedSequence(models.Model):
    user = models.ForeignKey(User, on_delete=models.CASCADE)
    sequence = models.TextField()
    created = models.DateTimeField(auto_now_add=True)

    def __str__(self):
        return f"{self.user.username} - {self.created}"
